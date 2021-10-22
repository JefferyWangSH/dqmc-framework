#include "Hubbard.h"
#include "CheckerBoard.h"
#include "SvdStack.h"
#include "BMatrixMult.hpp"
#include "StableGreens.hpp"
#include "Random.hpp"

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


Model::Hubbard::Hubbard(int ll, int lt, double beta, double t, double Uint, double mu, int nwrap, bool is_checkerboard)
{
    this->ll = ll;
    this->ls = ll * ll;
    this->lt = lt;

    this->beta = beta;
    this->dtau = beta / lt;
    this->Uint = Uint;
    this->alpha = acosh(exp(0.5 * dtau * abs(Uint)));
    this->u_is_attractive = (Uint < 0.0);
    this->is_checkerboard = is_checkerboard;

    this->t = t;
    this->mu = mu;

    this->nwrap = nwrap;
    this->current_tau = 0;

    // resize matrices and svdStacks
    this->s.resize(ls, lt);
    this->green_tt_up.resize(ls, ls);
    this->green_tt_dn.resize(ls, ls);
    this->green_t0_up.resize(ls, ls);
    this->green_t0_dn.resize(ls, ls);
    this->green_0t_up.resize(ls, ls);
    this->green_0t_dn.resize(ls, ls);

    this->stackLeftU = new SvdStack(ls, lt);
    this->stackLeftD = new SvdStack(ls, lt);
    this->stackRightU = new SvdStack(ls, lt);
    this->stackRightD = new SvdStack(ls, lt);

    this->vec_green_tt_up.reserve(lt);
    this->vec_green_tt_dn.reserve(lt);
    this->vec_green_t0_up.reserve(lt);
    this->vec_green_t0_dn.reserve(lt);
    this->vec_green_0t_up.reserve(lt);
    this->vec_green_0t_dn.reserve(lt);

    for (int l = 0; l < this->lt; ++l) {
        this->vec_green_tt_up.emplace_back(ls, ls);
        this->vec_green_tt_dn.emplace_back(ls, ls);
        this->vec_green_t0_up.emplace_back(ls, ls);
        this->vec_green_t0_dn.emplace_back(ls, ls);
        this->vec_green_0t_up.emplace_back(ls, ls);
        this->vec_green_0t_dn.emplace_back(ls, ls);
    }

    // set field config to random
    this->init_field_to_random();

    // initialize checkerboard
    this->checkerboard = new CheckerBoard::CheckerBoard();
    this->checkerboard->init_from_model(*this);
    this->is_checkerboard = this->checkerboard->is_checker_board();

    // initialize udv stacks for sweep use, stabilize every nwrap slices
    this->init_stacks(this->nwrap);

    // determine sign of current configuration
    this->config_sign = (this->green_tt_up.determinant() * this->green_tt_dn.determinant() >= 0)? +1.0 : -1.0;
}

Model::Hubbard::~Hubbard() {
    if (this->stackLeftU) {
        delete this->stackLeftU;
        this->stackLeftU = nullptr;
    }
    if (this->stackLeftD) {
        delete this->stackLeftD;
        this->stackLeftD = nullptr;
    }
    if (this->stackRightU) {
        delete this->stackRightU;
        this->stackRightU = nullptr;
    }
    if (this->stackRightD) {
        delete this->stackRightD;
        this->stackRightD = nullptr;
    }
    if (this->checkerboard) {
        delete this->checkerboard;
        this->checkerboard = nullptr;
    }
}

void Model::Hubbard::init_field_to_random() {
    // set field configuration to random
    assert( this->s.rows() == this->ls );
    assert( this->s.cols() == this->lt );

    std::bernoulli_distribution dist(0.5);
    for(int i = 0; i < this->ls; ++i) {
        for(int l = 0; l < this->lt; ++l) {
            this->s(i, l) = dist(Random::Engine)? +1.0:-1.0;
        }
    }
}

void Model::Hubbard::metropolis_update(int l) {
    /*
     * Update the aux boson field 's' at space-time position (i,l) for all i with Metropolis
     * probability, and - if the update is accepted - perform a
     * in-place update of the green's function.
     * Record the updated green's function at the life end of function.
     */
    assert( this->current_tau == l );
    assert( l >= 0 && l <= this->lt);

    const int tau = (l == 0)? this->lt-1 : l-1;
    for (int i = 0; i < this->ls; ++i) {
        // radio of flipping aux field s(i,l)
        double p;
        if (!this->u_is_attractive) {
            p = (1 + (1 - green_tt_up(i, i)) * (exp(-2 * alpha * s(i, tau)) - 1))
                       * (1 + (1 - green_tt_dn(i, i)) * (exp(+2 * alpha * s(i, tau)) - 1));
        }
        else {
            p = exp(2 * alpha * s(i, tau))  * (1 + (1 - green_tt_up(i, i)) * (exp(-2 * alpha * s(i, tau)) - 1))
                       * (1 + (1 - green_tt_dn(i, i)) * (exp(-2 * alpha * s(i, tau)) - 1));
        }

        if(std::bernoulli_distribution(std::min(1.0, std::abs(p)))(Random::Engine)) {
            /** reference:
             *  Quantum Monte Carlo Methods (Algorithms for Lattice Models) Determinant method
             *  Here we use the sparseness of matrix \delta */
            // update greens function (which is wrapped such that the update is at time slice 0 of g)
            // with a number of arithmetic operations proportional to N^2
            double factorU = (exp(-2 * alpha * s(i, tau)) - 1) / (1 + (1 - green_tt_up(i, i)) * (exp(-2 * alpha * s(i, tau)) - 1));
            this->green_tt_up -= factorU * this->green_tt_up.col(i) * (Eigen::VectorXd::Unit(ls, i).transpose() - this->green_tt_up.row(i));

            double factorD = (!this->u_is_attractive)? (exp(+2 * alpha * s(i, tau)) - 1) / (1 + (1 - green_tt_dn(i, i)) * (exp(+2 * alpha * s(i, tau)) - 1)) : factorU;
            this->green_tt_dn -= factorD * this->green_tt_dn.col(i) * (Eigen::VectorXd::Unit(ls, i).transpose() - this->green_tt_dn.row(i));

            // flip aux field
            this->s(i, tau) = -this->s(i, tau);

            // keep track of sign problem
            this->config_sign = (p >= 0)? +this->config_sign : -this->config_sign;
        }
    }
}

void Model::Hubbard::wrap_north(int l) {
    /*
     * Propagate the green's function from the current time slice l
     * upward to the time slice l+1:
     * G(l+1) = B_{l+1} G(l) B_{l+1}^{-1}
     * for both spin-1/2 state. Change green's functions in place.
     */
    assert( l >= 0 && l <= this->lt );

    const int tau = (l == this->lt)? 1 : l + 1;
    this->mult_B_from_left(this->green_tt_up, tau, +1);
    this->mult_invB_from_right(this->green_tt_up, tau, +1);
    this->mult_B_from_left(this->green_tt_dn, tau, -1);
    this->mult_invB_from_right(this->green_tt_dn, tau, -1);
}

void Model::Hubbard::wrap_south(int l) {
    /*
     * Propagate the green's function from the current time slice l
     * downward to the time slice l-1:
     * G(l-1) = B_{l}^{-1} G(l) B_{l}
     * for both spin-1/2 state. Change green's functions in place.
     */
    assert( l >= 0 && l <= this->lt );

    const int tau = (l == 0)? this->lt : l;
    this->mult_B_from_right(this->green_tt_up, tau, +1);
    this->mult_invB_from_left(this->green_tt_up, tau, +1);
    this->mult_B_from_right(this->green_tt_dn, tau, -1);
    this->mult_invB_from_left(this->green_tt_dn, tau, -1);
}

void Model::Hubbard::init_stacks(int stable_pace) {
    /*
     *  initialize udv stacks for sweep use
     *  sweep process will start from 0 to beta, so we initialize stackRight here.
     *  stabilize the process every stable_pace steps
     */
    assert( this->stackLeftU->empty() && this->stackLeftD->empty() );
    assert( this->stackRightU->empty() && this->stackRightD->empty() );

    Eigen::MatrixXd tmpU = Eigen::MatrixXd::Identity(ls, ls);
    Eigen::MatrixXd tmpD = Eigen::MatrixXd::Identity(ls, ls);

    // initial udv stacks for sweeping use
    for (int l = this->lt; l >= 1; --l) {
        this->mult_transB_from_left(tmpU, l, +1);
        this->mult_transB_from_left(tmpD, l, -1);

        // stabilize every stable_pace steps with svd decomposition
        if ((l - 1) % stable_pace == 0) {
            this->stackRightU->push(tmpU);
            this->stackRightD->push(tmpD);
            tmpU = Eigen::MatrixXd::Identity(ls, ls);
            tmpD = Eigen::MatrixXd::Identity(ls, ls);
        }
    }

    // initialize green function at l = 0
    compute_Green_eqtime(this->stackLeftU, this->stackRightU, this->green_tt_up);
    compute_Green_eqtime(this->stackLeftD, this->stackRightD, this->green_tt_dn);
}

void Model::Hubbard::sweep_0_to_beta(int stable_pace) {
    /*
     *  Update the space-time lattice of aux fields.
     *  For l = 1,2...,lt  flip fields and propagate green's functions
     *  Stabilize every `stable_pace` time slices
     */
    this->current_tau++;

    int nlen = (this->lt % stable_pace == 0)? this->lt/stable_pace : this->lt/stable_pace +1;
    assert( this->current_tau == 1 );
    assert( this->stackLeftU->empty() && this->stackLeftD->empty() );
    assert( this->stackRightU->len == nlen && this->stackRightD->len == nlen );

    // temporary matrices
    Eigen::MatrixXd tmpU = Eigen::MatrixXd::Identity(ls, ls);
    Eigen::MatrixXd tmpD = Eigen::MatrixXd::Identity(ls, ls);

    // sweep up from 0 to beta
    for (int l = 1; l <= this->lt; ++l) {
        // wrap green function to current time slice l
        this->wrap_north(l-1);

        // update aux field and record new greens
        this->metropolis_update(l);
        this->vec_green_tt_up[l-1] = this->green_tt_up;
        this->vec_green_tt_dn[l-1] = this->green_tt_dn;

        this->mult_B_from_left(tmpU, l, +1);
        this->mult_B_from_left(tmpD, l, -1);

        if (l % stable_pace == 0 || l == this->lt) {
            // wrap greens function
            this->stackRightU->pop();
            this->stackRightD->pop();
            this->stackLeftU->push(tmpU);
            this->stackLeftD->push(tmpD);

            Eigen::MatrixXd tmp_green_tt_up = Eigen::MatrixXd::Zero(ls, ls);
            Eigen::MatrixXd tmp_green_tt_dn = Eigen::MatrixXd::Zero(ls, ls);
            double tmp_wrap_error_tt_up = 0.0;
            double tmp_wrap_error_tt_dn = 0.0;

            // compute fresh greens every `stable_pace` steps: g = (1 + stackLeft * stackRight^T)^-1
            // stackLeft = B(l-1) *...* B(0)
            // stackRight = B(l)^T *...* B(L-1)^T
            compute_Green_eqtime(this->stackLeftU, this->stackRightU, tmp_green_tt_up);
            compute_Green_eqtime(this->stackLeftD, this->stackRightD, tmp_green_tt_dn);

            // calculate wrap error
            matrix_compare_error(tmp_green_tt_up, this->green_tt_up, tmp_wrap_error_tt_up);
            matrix_compare_error(tmp_green_tt_dn, this->green_tt_dn, tmp_wrap_error_tt_dn);
            this->max_wrap_error_equal = std::max(this->max_wrap_error_equal, std::max(tmp_wrap_error_tt_up, tmp_wrap_error_tt_dn));

            this->green_tt_up = tmp_green_tt_up;
            this->green_tt_dn = tmp_green_tt_dn;

            this->vec_green_tt_up[l-1] = this->green_tt_up;
            this->vec_green_tt_dn[l-1] = this->green_tt_dn;

            tmpU = Eigen::MatrixXd::Identity(ls, ls);
            tmpD = Eigen::MatrixXd::Identity(ls, ls);
        }

        // in the end stop at l = lt + 1
        this->current_tau++;
    }

    // end with fresh green functions
    this->vec_green_tt_up[lt-1] = this->green_tt_up;
    this->vec_green_tt_dn[lt-1] = this->green_tt_dn;
}

void Model::Hubbard::sweep_beta_to_0(int stable_pace) {
    /*
     *  Update the space-time lattice of aux fields.
     *  For l=lt,lt-1,...,1  flip fields and propagate green's functions
     *  Stabilize every `stable_pace` time slices
     */
    this->current_tau--;

    int nlen = (this->lt % stable_pace == 0)? this->lt/stable_pace : this->lt/stable_pace + 1;
    assert( this->current_tau == this->lt );
    assert( this->stackRightU->empty() && this->stackRightD->empty() );
    assert( this->stackLeftU->len == nlen && this->stackLeftD->len == nlen );

    // temporary matrices
    Eigen::MatrixXd tmpU = Eigen::MatrixXd::Identity(ls, ls);
    Eigen::MatrixXd tmpD = Eigen::MatrixXd::Identity(ls, ls);

    // sweep down from beta to 0
    for (int l = this->lt; l >= 1; --l) {
        if (l % stable_pace == 0 && l != this->lt) {
            // update udv stacks
            this->stackLeftU->pop();
            this->stackLeftD->pop();
            this->stackRightU->push(tmpU);
            this->stackRightD->push(tmpD);

            Eigen::MatrixXd tmp_green_tt_up = Eigen::MatrixXd::Zero(ls, ls);
            Eigen::MatrixXd tmp_green_tt_dn = Eigen::MatrixXd::Zero(ls, ls);
            double tmp_wrap_error_tt_up = 0.0;
            double tmp_wrap_error_tt_dn = 0.0;

            compute_Green_eqtime(this->stackLeftU, this->stackRightU, tmp_green_tt_up);
            compute_Green_eqtime(this->stackLeftD, this->stackRightD, tmp_green_tt_dn);

            // calculate wrap error
            matrix_compare_error(tmp_green_tt_up, this->green_tt_up, tmp_wrap_error_tt_up);
            matrix_compare_error(tmp_green_tt_dn, this->green_tt_dn, tmp_wrap_error_tt_dn);
            this->max_wrap_error_equal = std::max(this->max_wrap_error_equal, std::max(tmp_wrap_error_tt_up, tmp_wrap_error_tt_dn));

            this->green_tt_up = tmp_green_tt_up;
            this->green_tt_dn = tmp_green_tt_dn;

            tmpU = Eigen::MatrixXd::Identity(ls, ls);
            tmpD = Eigen::MatrixXd::Identity(ls, ls);
        }

        // update aux field and record new greens
        this->metropolis_update(l);
        this->vec_green_tt_up[l-1] = this->green_tt_up;
        this->vec_green_tt_dn[l-1] = this->green_tt_dn;

        this->mult_transB_from_left(tmpU, l, +1);
        this->mult_transB_from_left(tmpD, l, -1);

        this->wrap_south(l);

        this->current_tau--;
    }

    // at l = 0
    this->stackLeftU->pop();
    this->stackLeftD->pop();
    this->stackRightU->push(tmpU);
    this->stackRightD->push(tmpD);

    compute_Green_eqtime(this->stackLeftU, this->stackRightU, this->green_tt_up);
    compute_Green_eqtime(this->stackLeftD, this->stackRightD, this->green_tt_dn);

    // end with fresh green functions
    this->vec_green_tt_up[lt-1] = this->green_tt_up;
    this->vec_green_tt_dn[lt-1] = this->green_tt_dn;
}

void Model::Hubbard::sweep_0_to_beta_displaced(int stable_pace) {
    /*
     *  Calculate time-displaced green's function, while the aux field remains unchanged.
     *  For l = 1,2...,lt, recompute SvdStacks every `stable_pace` time slices.
     *  Data is stored in vec_green_t0/0t_up/dn.
     *
     *  Cautious that equal-time green's functions are also re-calculated according to the current configurations of aux field.
     *  Data is stored in vec_green_tt_up/dn
     */
    this->current_tau++;

    int nlen = (this->lt % stable_pace == 0)? this->lt/stable_pace : this->lt/stable_pace + 1;
    assert( this->current_tau == 1 );
    assert( this->stackLeftU->empty() && this->stackLeftD->empty() );
    assert( this->stackRightU->len == nlen && this->stackRightD->len == nlen );

    // initialize: at l = 0, gt0 = g00, g0t = g00 - 1
    this->green_t0_up = this->green_tt_up;
    this->green_t0_dn = this->green_tt_dn;
    this->green_0t_up = this->green_tt_up - Eigen::MatrixXd::Identity(ls, ls);
    this->green_0t_dn = this->green_tt_dn - Eigen::MatrixXd::Identity(ls, ls);

    // temporary matrices
    Eigen::MatrixXd tmpU = Eigen::MatrixXd::Identity(ls, ls);
    Eigen::MatrixXd tmpD = Eigen::MatrixXd::Identity(ls, ls);

    // sweep up from 0 to beta
    for (int l = 1; l <= this->lt; ++l) {
        // wrap equal time green function to current time slice l
        this->wrap_north(l-1);
        this->vec_green_tt_up[l-1] = this->green_tt_up;
        this->vec_green_tt_dn[l-1] = this->green_tt_dn;
        
        // calculate and record time-displaced green functions at different time slices
        this->mult_B_from_left(this->green_t0_up, l, +1);
        this->mult_B_from_left(this->green_t0_dn, l, -1);
        this->vec_green_t0_up[l-1] = this->green_t0_up;
        this->vec_green_t0_dn[l-1] = this->green_t0_dn;

        this->mult_invB_from_right(this->green_0t_up, l, +1);
        this->mult_invB_from_right(this->green_0t_dn, l, -1);
        this->vec_green_0t_up[l-1] = this->green_0t_up;
        this->vec_green_0t_dn[l-1] = this->green_0t_dn;

        this->mult_B_from_left(tmpU, l, +1);
        this->mult_B_from_left(tmpD, l, -1);

        if (l % stable_pace == 0 || l == this->lt) {
            // wrap greens function
            this->stackRightU->pop();
            this->stackRightD->pop();
            this->stackLeftU->push(tmpU);
            this->stackLeftD->push(tmpD);

            Eigen::MatrixXd tmp_green_t0_up = Eigen::MatrixXd::Zero(ls, ls);
            Eigen::MatrixXd tmp_green_t0_dn = Eigen::MatrixXd::Zero(ls, ls);
            Eigen::MatrixXd tmp_green_0t_up = Eigen::MatrixXd::Zero(ls, ls);
            Eigen::MatrixXd tmp_green_0t_dn = Eigen::MatrixXd::Zero(ls, ls);
            double tmp_wrap_error_t0_up = 0.0;
            double tmp_wrap_error_t0_dn = 0.0;
            double tmp_wrap_error_0t_up = 0.0;
            double tmp_wrap_error_0t_dn = 0.0;

            // compute fresh greens every stable_pace steps
            // stackLeft = B(l-1) *...* B(0)
            // stackRight = B(l)^T *...* B(L-1)^T
            // equal time green's function are re-evaluated for current field configurations
            compute_Green_eqtime(this->stackLeftU, this->stackRightU, this->green_tt_up);
            compute_Green_eqtime(this->stackLeftD, this->stackRightD, this->green_tt_dn);
            compute_Green_displaced(this->stackLeftU, this->stackRightU, tmp_green_t0_up, tmp_green_0t_up);
            compute_Green_displaced(this->stackLeftD, this->stackRightD, tmp_green_t0_dn, tmp_green_0t_dn);

            // calculate wrap error
            matrix_compare_error(tmp_green_t0_up, this->green_t0_up, tmp_wrap_error_t0_up);
            matrix_compare_error(tmp_green_t0_dn, this->green_t0_dn, tmp_wrap_error_t0_dn);
            this->max_wrap_error_displaced = std::max(this->max_wrap_error_displaced, std::max(tmp_wrap_error_t0_up, tmp_wrap_error_t0_dn));

            matrix_compare_error(tmp_green_0t_up, this->green_0t_up, tmp_wrap_error_0t_up);
            matrix_compare_error(tmp_green_t0_dn, this->green_t0_dn, tmp_wrap_error_t0_dn);
            this->max_wrap_error_displaced = std::max(this->max_wrap_error_displaced, std::max(tmp_wrap_error_0t_up, tmp_wrap_error_0t_dn));

            this->green_t0_up = tmp_green_t0_up;
            this->green_t0_dn = tmp_green_t0_dn;
            this->green_0t_up = tmp_green_0t_up;
            this->green_0t_dn = tmp_green_0t_dn;

            this->vec_green_tt_up[l-1] = this->green_tt_up;
            this->vec_green_tt_dn[l-1] = this->green_tt_dn;
            this->vec_green_t0_up[l-1] = this->green_t0_up;
            this->vec_green_t0_dn[l-1] = this->green_t0_dn;
            this->vec_green_0t_up[l-1] = this->green_0t_up;
            this->vec_green_0t_dn[l-1] = this->green_0t_dn;

            tmpU = Eigen::MatrixXd::Identity(ls, ls);
            tmpD = Eigen::MatrixXd::Identity(ls, ls);
        }

        // in the end stop at l = lt + 1
        this->current_tau++;
    }
}
