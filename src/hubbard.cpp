#include "hubbard.h"
#include "checker_board.h"
#include "svd_stack.h"
#include "b_matrix_mult.hpp"
#include "stable_greens.hpp"
#include "random.h"

#include <iostream>
#include <cassert>

/** TODO:
  *   1. cyclic sweeping (done)
  *   2. chemical potential mu (done)
  *   3. stable multiplication of ill-conditioned matrices (done)
  *   4. check wrap error of green's function (done)
  *   6. detQMC class (done)
  *   7. checkerboard decomposition (done)
  *   8. time-displaced dynamical measurements: green_t0/0t (done)
  *   9. attractive interaction U < 0 (done)
  *   10. reweighing for doped case (done)
  *   11. ...
  */


void Model::Hubbard::set_model_params(int ll, int lt, double beta, double t, double u_int, double mu) {
    this->ll = ll;
    this->ls = ll * ll;
    this->lt = lt;

    this->beta = beta;
    this->dtau = beta / lt;
    this->u_int = u_int;
    this->alpha = acosh(exp(0.5 * dtau * abs(u_int)));
    this->is_attractive_u = (u_int < 0.0);
    this->t = t;
    this->mu = mu;

    this->current_tau = 0;
}

void Model::Hubbard::set_bool_params(bool is_eqtime_measure, bool is_dynamic_measure, bool is_checkerboard) {
    this->is_eqtime_measure = is_eqtime_measure;
    this->is_dynamic_measure = is_dynamic_measure;
    this->is_checkerboard = is_checkerboard;
}

void Model::Hubbard::set_stabilization_pace(int nwrap) {
    this->nwrap = nwrap;
}

void Model::Hubbard::deallocate() {
    if (this->s) { s.reset(); }
    if (this->green_tt_up) { this->green_tt_up.reset(); }
    if (this->green_tt_dn) { this->green_tt_dn.reset(); }
    if (this->green_t0_up) { this->green_t0_up.reset(); }
    if (this->green_t0_dn) { this->green_t0_dn.reset(); }
    if (this->green_0t_up) { this->green_0t_up.reset(); }
    if (this->green_0t_dn) { this->green_0t_dn.reset(); }
    if (this->vec_green_tt_up) { this->vec_green_tt_up.reset(); }
    if (this->vec_green_tt_dn) { this->vec_green_tt_dn.reset(); }
    if (this->vec_green_t0_up) { this->vec_green_t0_up.reset(); }
    if (this->vec_green_t0_dn) { this->vec_green_t0_dn.reset(); }
    if (this->vec_green_0t_up) { this->vec_green_0t_up.reset(); }
    if (this->vec_green_0t_dn) { this->vec_green_0t_dn.reset(); }
    if (this->stack_left_up) { this->stack_left_up.reset(); }
    if (this->stack_left_dn) { this->stack_left_dn.reset(); }
    if (this->stack_right_up) { this->stack_right_up.reset(); }
    if (this->stack_right_dn) { this->stack_right_dn.reset(); }
    if (this->checkerboard) { this->checkerboard.reset(); }
}

void Model::Hubbard::allocate() {
    // release memory if allocated before
    this->deallocate();

    // resize matrices and svdStacks
    this->s = std::make_unique<Eigen::MatrixXd>(this->ls, this->lt);
    this->green_tt_up = std::make_unique<Eigen::MatrixXd>(this->ls, this->ls);
    this->green_tt_dn = std::make_unique<Eigen::MatrixXd>(this->ls, this->ls);
    if (this->is_eqtime_measure) {
        this->vec_green_tt_up = std::make_unique<std::vector<Eigen::MatrixXd>>(this->lt, Eigen::MatrixXd(this->ls, this->ls));
        this->vec_green_tt_dn = std::make_unique<std::vector<Eigen::MatrixXd>>(this->lt, Eigen::MatrixXd(this->ls, this->ls));
        this->vec_config_sign = std::make_unique<std::vector<double>>(this->lt, 0.0);
    }

    if (this->is_dynamic_measure) {
        this->green_t0_up = std::make_unique<Eigen::MatrixXd>(this->ls, this->ls);
        this->green_t0_dn = std::make_unique<Eigen::MatrixXd>(this->ls, this->ls);
        this->green_0t_up = std::make_unique<Eigen::MatrixXd>(this->ls, this->ls);
        this->green_0t_dn = std::make_unique<Eigen::MatrixXd>(this->ls, this->ls);
        this->vec_green_t0_up = std::make_unique<std::vector<Eigen::MatrixXd>>(this->lt, Eigen::MatrixXd(this->ls, this->ls));
        this->vec_green_t0_dn = std::make_unique<std::vector<Eigen::MatrixXd>>(this->lt, Eigen::MatrixXd(this->ls, this->ls));
        this->vec_green_0t_up = std::make_unique<std::vector<Eigen::MatrixXd>>(this->lt, Eigen::MatrixXd(this->ls, this->ls));
        this->vec_green_0t_dn = std::make_unique<std::vector<Eigen::MatrixXd>>(this->lt, Eigen::MatrixXd(this->ls, this->ls));
    }

    this->stack_left_up = std::make_unique<SvdStack>(this->ls, this->lt);
    this->stack_left_dn = std::make_unique<SvdStack>(this->ls, this->lt);
    this->stack_right_up = std::make_unique<SvdStack>(this->ls, this->lt);
    this->stack_right_dn = std::make_unique<SvdStack>(this->ls, this->lt);

    this->checkerboard = std::make_unique<CheckerBoard::CheckerBoard>();
}

void Model::Hubbard::init_field_to_random() {
    // set field configuration to random
    assert( this->s->rows() == this->ls );
    assert( this->s->cols() == this->lt );

    std::bernoulli_distribution dist(0.5);
    for(int i = 0; i < this->ls; ++i) {
        for(int l = 0; l < this->lt; ++l) {
            (*this->s)(i, l) = dist(Random::Engine)? +1.0 : -1.0;
        }
    }
}

void Model::Hubbard::initial() {
    // allocate memory
    this->allocate();

    // initialize checkerboard
    this->checkerboard->init_from_model(*this);
    this->is_checkerboard = this->checkerboard->is_checker_board();

    // set field config to random
    this->init_field_to_random();

    // initialize udv stacks for sweep use, stabilize every nwrap slices
    this->init_stacks();

    // determine sign of current configuration
    this->config_sign = (this->green_tt_up->determinant() * this->green_tt_dn->determinant() >= 0)? +1.0 : -1.0;
}

void Model::Hubbard::metropolis_update(int l) {
    /*
     * Update the aux boson field 's' at space-time position (i,l) for all i with Metropolis
     * probability, and - if the update is accepted - perform a
     * in-place update of the green's function.
     * Record the updated green's function at the life end of function.
     */
    assert( this->current_tau == l );
    assert( l >= 0 && l <= this->lt );

    const int tau = (l == 0)? this->lt-1 : l-1;
    for (int i = 0; i < this->ls; ++i) {
        // radio of flipping aux field s(i,l)
        double p;
        if (!this->is_attractive_u) {
            p =   (1 + (1 - (*this->green_tt_up)(i, i)) * (exp(-2 * this->alpha * (*this->s)(i, tau)) - 1))
                * (1 + (1 - (*this->green_tt_dn)(i, i)) * (exp(+2 * this->alpha * (*this->s)(i, tau)) - 1));
        }
        else {
            p = exp(2 * this->alpha * (*this->s)(i, tau))
                * (1 + (1 - (*this->green_tt_up)(i, i)) * (exp(-2 * this->alpha * (*this->s)(i, tau)) - 1))
                * (1 + (1 - (*this->green_tt_dn)(i, i)) * (exp(-2 * this->alpha * (*this->s)(i, tau)) - 1));
        }

        if(std::bernoulli_distribution(std::min(1.0, std::abs(p)))(Random::Engine)) {
            /** reference:
             *  Quantum Monte Carlo Methods (Algorithms for Lattice Models) Determinant method
             *  Here we use the sparseness of matrix \delta */
            // update greens function (which is wrapped such that the update is at time slice 0 of g)
            // with a number of arithmetic operations proportional to N^2
            double factor_up = (exp(-2 * this->alpha * (*this->s)(i, tau)) - 1)
                             / (1 + (1 - (*this->green_tt_up)(i, i)) * (exp(-2 * this->alpha * (*this->s)(i, tau)) - 1));
            (*this->green_tt_up) -= factor_up * (*this->green_tt_up).col(i) 
                                              * (Eigen::VectorXd::Unit(this->ls, i).transpose() - (*this->green_tt_up).row(i));
            double factor_dn = (!this->is_attractive_u)? 
                               (exp(+2 * this->alpha * (*this->s)(i, tau)) - 1) / (1 + (1 - (*this->green_tt_dn)(i, i))
                             * (exp(+2 * this->alpha * (*this->s)(i, tau)) - 1)) : factor_up;
            (*this->green_tt_dn) -= factor_dn * this->green_tt_dn->col(i)
                                              * (Eigen::VectorXd::Unit(this->ls, i).transpose() - this->green_tt_dn->row(i));

            // flip aux field
            (*this->s)(i, tau) = -(*this->s)(i, tau);

            // keep track of sign problem
            this->config_sign = (p >= 0)? +this->config_sign : -this->config_sign;
        }
    }
}

void Model::Hubbard::wrap_0_to_beta(int l) {
    /*
     * Propagate the green's function from the current time slice l
     * upward to the time slice l+1:
     * G(l+1) = B_{l+1} G(l) B_{l+1}^{-1}
     * for both spin-1/2 state. Change green's functions in place.
     */
    assert( l >= 0 && l <= this->lt );

    const int tau = (l == this->lt)? 1 : l + 1;
    this->mult_B_from_left(*this->green_tt_up, tau, +1);
    this->mult_invB_from_right(*this->green_tt_up, tau, +1);
    this->mult_B_from_left(*this->green_tt_dn, tau, -1);
    this->mult_invB_from_right(*this->green_tt_dn, tau, -1);
}

void Model::Hubbard::wrap_beta_to_0(int l) {
    /*
     * Propagate the green's function from the current time slice l
     * downward to the time slice l-1:
     * G(l-1) = B_{l}^{-1} G(l) B_{l}
     * for both spin-1/2 state. Change green's functions in place.
     */
    assert( l >= 0 && l <= this->lt );

    const int tau = (l == 0)? this->lt : l;
    this->mult_B_from_right(*this->green_tt_up, tau, +1);
    this->mult_invB_from_left(*this->green_tt_up, tau, +1);
    this->mult_B_from_right(*this->green_tt_dn, tau, -1);
    this->mult_invB_from_left(*this->green_tt_dn, tau, -1);
}

void Model::Hubbard::init_stacks() {
    /*
     *  initialize udv stacks for sweep use
     *  sweep process will start from 0 to beta, so we initialize stack_right here.
     *  stabilize the process every nwrap steps
     */
    assert( this->stack_left_up->empty() && this->stack_left_dn->empty() );
    assert( this->stack_right_up->empty() && this->stack_right_dn->empty() );

    Eigen::MatrixXd tmp_up = Eigen::MatrixXd::Identity(this->ls, this->ls);
    Eigen::MatrixXd tmp_dn = Eigen::MatrixXd::Identity(this->ls, this->ls);

    // initial udv stacks for sweeping use
    for (int l = this->lt; l >= 1; --l) {
        this->mult_transB_from_left(tmp_up, l, +1);
        this->mult_transB_from_left(tmp_dn, l, -1);

        // stabilize every nwrap steps with svd decomposition
        if ((l - 1) % this->nwrap == 0) {
            this->stack_right_up->push(tmp_up);
            this->stack_right_dn->push(tmp_dn);
            tmp_up = Eigen::MatrixXd::Identity(this->ls, this->ls);
            tmp_dn = Eigen::MatrixXd::Identity(this->ls, this->ls);
        }
    }

    // initialize green function at l = 0
    StableGreens::compute_greens_eqtime(this->stack_left_up.get(), this->stack_right_up.get(), *this->green_tt_up);
    StableGreens::compute_greens_eqtime(this->stack_left_dn.get(), this->stack_right_dn.get(), *this->green_tt_dn);
}

void Model::Hubbard::sweep_0_to_beta() {
    /*
     *  Update the space-time lattice of aux fields.
     *  For l = 1,2...,lt  flip fields and propagate green's functions
     *  Stabilize every nwrap time slices
     */
    this->current_tau++;

    int nlen = (this->lt % this->nwrap == 0)? this->lt/this->nwrap : this->lt/this->nwrap +1;
    assert( this->current_tau == 1 );
    assert( this->stack_left_up->empty() && this->stack_left_dn->empty() );
    assert( this->stack_right_up->length() == nlen && this->stack_right_dn->length() == nlen );

    // temporary matrices
    Eigen::MatrixXd tmp_up = Eigen::MatrixXd::Identity(this->ls, this->ls);
    Eigen::MatrixXd tmp_dn = Eigen::MatrixXd::Identity(this->ls, this->ls);

    // sweep up from 0 to beta
    for (int l = 1; l <= this->lt; ++l) {
        // wrap green function to current time slice l
        this->wrap_0_to_beta(l-1);

        // update aux field and record new greens
        this->metropolis_update(l);
        if (this->is_eqtime_measure) {
            (*this->vec_green_tt_up)[l-1] = *this->green_tt_up;
            (*this->vec_green_tt_dn)[l-1] = *this->green_tt_dn;
            (*this->vec_config_sign)[l-1] = this->config_sign;
        }

        this->mult_B_from_left(tmp_up, l, +1);
        this->mult_B_from_left(tmp_dn, l, -1);

        if (l % this->nwrap == 0 || l == this->lt) {
            // wrap greens function
            this->stack_right_up->pop();
            this->stack_right_dn->pop();
            this->stack_left_up->push(tmp_up);
            this->stack_left_dn->push(tmp_dn);

            Eigen::MatrixXd tmp_green_tt_up = Eigen::MatrixXd::Zero(this->ls, this->ls);
            Eigen::MatrixXd tmp_green_tt_dn = Eigen::MatrixXd::Zero(this->ls, this->ls);
            double tmp_wrap_error_tt_up = 0.0;
            double tmp_wrap_error_tt_dn = 0.0;

            // compute fresh greens every nwrap steps: g = (1 + stack_left * stack_right^T)^-1
            // stack_left = B(l-1) *...* B(0)
            // stack_right = B(l)^T *...* B(L-1)^T
            StableGreens::compute_greens_eqtime(this->stack_left_up.get(), this->stack_right_up.get(), tmp_green_tt_up);
            StableGreens::compute_greens_eqtime(this->stack_left_dn.get(), this->stack_right_dn.get(), tmp_green_tt_dn);

            // calculate wrap error
            StableGreens::matrix_compare_error(tmp_green_tt_up, *this->green_tt_up, tmp_wrap_error_tt_up);
            StableGreens::matrix_compare_error(tmp_green_tt_dn, *this->green_tt_dn, tmp_wrap_error_tt_dn);
            this->max_wrap_error_equal = std::max(this->max_wrap_error_equal, std::max(tmp_wrap_error_tt_up, tmp_wrap_error_tt_dn));

            *this->green_tt_up = tmp_green_tt_up;
            *this->green_tt_dn = tmp_green_tt_dn;

            if (this->is_eqtime_measure) {
                (*this->vec_green_tt_up)[l-1] = *this->green_tt_up;
                (*this->vec_green_tt_dn)[l-1] = *this->green_tt_dn;
            }

            tmp_up = Eigen::MatrixXd::Identity(this->ls, this->ls);
            tmp_dn = Eigen::MatrixXd::Identity(this->ls, this->ls);
        }

        // finally stop at l = lt + 1
        this->current_tau++;
    }

    // end with fresh green functions
    if (this->is_eqtime_measure) {
        (*this->vec_green_tt_up)[lt-1] = *this->green_tt_up;
        (*this->vec_green_tt_dn)[lt-1] = *this->green_tt_dn;
    }
}

void Model::Hubbard::sweep_beta_to_0() {
    /*
     *  Update the space-time lattice of aux fields.
     *  For l=lt,lt-1,...,1  flip fields and propagate green's functions
     *  Stabilize every nwrap time slices
     */
    this->current_tau--;

    int nlen = (this->lt % this->nwrap == 0)? this->lt/this->nwrap : this->lt/this->nwrap + 1;
    assert( this->current_tau == this->lt );
    assert( this->stack_right_up->empty() && this->stack_right_dn->empty() );
    assert( this->stack_left_up->length() == nlen && this->stack_left_dn->length() == nlen );

    // temporary matrices
    Eigen::MatrixXd tmp_up = Eigen::MatrixXd::Identity(this->ls, this->ls);
    Eigen::MatrixXd tmp_dn = Eigen::MatrixXd::Identity(this->ls, this->ls);

    // sweep down from beta to 0
    for (int l = this->lt; l >= 1; --l) {
        if (l % this->nwrap == 0 && l != this->lt) {
            // update udv stacks
            this->stack_left_up->pop();
            this->stack_left_dn->pop();
            this->stack_right_up->push(tmp_up);
            this->stack_right_dn->push(tmp_dn);

            Eigen::MatrixXd tmp_green_tt_up = Eigen::MatrixXd::Zero(this->ls, this->ls);
            Eigen::MatrixXd tmp_green_tt_dn = Eigen::MatrixXd::Zero(this->ls, this->ls);
            double tmp_wrap_error_tt_up = 0.0;
            double tmp_wrap_error_tt_dn = 0.0;

            StableGreens::compute_greens_eqtime(this->stack_left_up.get(), this->stack_right_up.get(), tmp_green_tt_up);
            StableGreens::compute_greens_eqtime(this->stack_left_dn.get(), this->stack_right_dn.get(), tmp_green_tt_dn);

            // calculate wrap error
            StableGreens::matrix_compare_error(tmp_green_tt_up, *this->green_tt_up, tmp_wrap_error_tt_up);
            StableGreens::matrix_compare_error(tmp_green_tt_dn, *this->green_tt_dn, tmp_wrap_error_tt_dn);
            this->max_wrap_error_equal = std::max(this->max_wrap_error_equal, std::max(tmp_wrap_error_tt_up, tmp_wrap_error_tt_dn));

            *this->green_tt_up = tmp_green_tt_up;
            *this->green_tt_dn = tmp_green_tt_dn;

            tmp_up = Eigen::MatrixXd::Identity(this->ls, this->ls);
            tmp_dn = Eigen::MatrixXd::Identity(this->ls, this->ls);
        }

        // update aux field and record new greens
        this->metropolis_update(l);
        if (this->is_eqtime_measure) {
            (*this->vec_green_tt_up)[l-1] = *this->green_tt_up;
            (*this->vec_green_tt_dn)[l-1] = *this->green_tt_dn;
            (*this->vec_config_sign)[l-1] = this->config_sign;
        }

        this->mult_transB_from_left(tmp_up, l, +1);
        this->mult_transB_from_left(tmp_dn, l, -1);

        this->wrap_beta_to_0(l);

        this->current_tau--;
    }

    // at l = 0
    this->stack_left_up->pop();
    this->stack_left_dn->pop();
    this->stack_right_up->push(tmp_up);
    this->stack_right_dn->push(tmp_dn);

    StableGreens::compute_greens_eqtime(this->stack_left_up.get(), this->stack_right_up.get(), *this->green_tt_up);
    StableGreens::compute_greens_eqtime(this->stack_left_dn.get(), this->stack_right_dn.get(), *this->green_tt_dn);

    // end with fresh green functions
    if (this->is_eqtime_measure) {
        (*this->vec_green_tt_up)[lt-1] = *this->green_tt_up;
        (*this->vec_green_tt_dn)[lt-1] = *this->green_tt_dn;
    }
}

void Model::Hubbard::sweep_0_to_beta_dynamic() {
    /*
     *  Calculate time-displaced (dynamical) green's function, while the aux field remains unchanged.
     *  For l = 1,2...,lt, recompute SvdStacks every nwrap time slices.
     *  Data is stored in vec_green_t0/0t_up/dn.
     *
     *  Cautious that equal-time green's functions are also re-calculated according to the current configurations of aux field.
     *  Data is stored in vec_green_tt_up/dn
     */
    if (this->is_dynamic_measure) {

        this->current_tau++;
        int nlen = (this->lt % this->nwrap == 0)? this->lt/this->nwrap : this->lt/this->nwrap + 1;
        assert( this->current_tau == 1 );
        assert( this->stack_left_up->empty() && this->stack_left_dn->empty() );
        assert( this->stack_right_up->length() == nlen && this->stack_right_dn->length() == nlen );

        // initialize: at l = 0, gt0 = g00, g0t = g00 - 1
        *this->green_t0_up = *this->green_tt_up;
        *this->green_t0_dn = *this->green_tt_dn;
        *this->green_0t_up = *this->green_tt_up - Eigen::MatrixXd::Identity(this->ls, this->ls);
        *this->green_0t_dn = *this->green_tt_dn - Eigen::MatrixXd::Identity(this->ls, this->ls);

        // temporary matrices
        Eigen::MatrixXd tmp_up = Eigen::MatrixXd::Identity(this->ls, this->ls);
        Eigen::MatrixXd tmp_dn = Eigen::MatrixXd::Identity(this->ls, this->ls);

        // sweep up from 0 to beta
        for (int l = 1; l <= this->lt; ++l) {
            // wrap equal time green function to current time slice l
            this->wrap_0_to_beta(l-1);
            (*this->vec_green_tt_up)[l-1] = *this->green_tt_up;
            (*this->vec_green_tt_dn)[l-1] = *this->green_tt_dn;
            
            // calculate and record time-displaced green functions at different time slices
            this->mult_B_from_left(*this->green_t0_up, l, +1);
            this->mult_B_from_left(*this->green_t0_dn, l, -1);
            (*this->vec_green_t0_up)[l-1] = *this->green_t0_up;
            (*this->vec_green_t0_dn)[l-1] = *this->green_t0_dn;

            this->mult_invB_from_right(*this->green_0t_up, l, +1);
            this->mult_invB_from_right(*this->green_0t_dn, l, -1);
            (*this->vec_green_0t_up)[l-1] = *this->green_0t_up;
            (*this->vec_green_0t_dn)[l-1] = *this->green_0t_dn;

            this->mult_B_from_left(tmp_up, l, +1);
            this->mult_B_from_left(tmp_dn, l, -1);

            if (l % this->nwrap == 0 || l == this->lt) {
                // wrap greens function
                this->stack_right_up->pop();
                this->stack_right_dn->pop();
                this->stack_left_up->push(tmp_up);
                this->stack_left_dn->push(tmp_dn);

                Eigen::MatrixXd tmp_green_t0_up = Eigen::MatrixXd::Zero(this->ls, this->ls);
                Eigen::MatrixXd tmp_green_t0_dn = Eigen::MatrixXd::Zero(this->ls, this->ls);
                Eigen::MatrixXd tmp_green_0t_up = Eigen::MatrixXd::Zero(this->ls, this->ls);
                Eigen::MatrixXd tmp_green_0t_dn = Eigen::MatrixXd::Zero(this->ls, this->ls);
                double tmp_wrap_error_t0_up = 0.0;
                double tmp_wrap_error_t0_dn = 0.0;
                double tmp_wrap_error_0t_up = 0.0;
                double tmp_wrap_error_0t_dn = 0.0;

                // compute fresh greens every nwrap steps
                // stack_left = B(l-1) *...* B(0)
                // stack_right = B(l)^T *...* B(L-1)^T
                // equal time green's function are re-evaluated for current field configurations
                StableGreens::compute_greens_eqtime(this->stack_left_up.get(), this->stack_right_up.get(), *this->green_tt_up);
                StableGreens::compute_greens_eqtime(this->stack_left_dn.get(), this->stack_right_dn.get(), *this->green_tt_dn);
                StableGreens::compute_greens_dynamic(this->stack_left_up.get(), this->stack_right_up.get(), tmp_green_t0_up, tmp_green_0t_up);
                StableGreens::compute_greens_dynamic(this->stack_left_dn.get(), this->stack_right_dn.get(), tmp_green_t0_dn, tmp_green_0t_dn);

                // calculate wrap error
                StableGreens::matrix_compare_error(tmp_green_t0_up, *this->green_t0_up, tmp_wrap_error_t0_up);
                StableGreens::matrix_compare_error(tmp_green_t0_dn, *this->green_t0_dn, tmp_wrap_error_t0_dn);
                this->max_wrap_error_dynamic = std::max(this->max_wrap_error_dynamic, std::max(tmp_wrap_error_t0_up, tmp_wrap_error_t0_dn));

                StableGreens::matrix_compare_error(tmp_green_0t_up, *this->green_0t_up, tmp_wrap_error_0t_up);
                StableGreens::matrix_compare_error(tmp_green_0t_dn, *this->green_0t_dn, tmp_wrap_error_0t_dn);
                this->max_wrap_error_dynamic = std::max(this->max_wrap_error_dynamic, std::max(tmp_wrap_error_0t_up, tmp_wrap_error_0t_dn));

                *this->green_t0_up = tmp_green_t0_up;
                *this->green_t0_dn = tmp_green_t0_dn;
                *this->green_0t_up = tmp_green_0t_up;
                *this->green_0t_dn = tmp_green_0t_dn;

                (*this->vec_green_tt_up)[l-1] = *this->green_tt_up;
                (*this->vec_green_tt_dn)[l-1] = *this->green_tt_dn;
                (*this->vec_green_t0_up)[l-1] = *this->green_t0_up;
                (*this->vec_green_t0_dn)[l-1] = *this->green_t0_dn;
                (*this->vec_green_0t_up)[l-1] = *this->green_0t_up;
                (*this->vec_green_0t_dn)[l-1] = *this->green_0t_dn;

                tmp_up = Eigen::MatrixXd::Identity(this->ls, this->ls);
                tmp_dn = Eigen::MatrixXd::Identity(this->ls, this->ls);
            }

            // finally stop at l = lt + 1
            this->current_tau++;
        }
    }
}
