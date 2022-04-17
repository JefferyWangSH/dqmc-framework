#include "dqmc_walker.h"
#include "svd_stack.h"
#include "lattice/lattice_base.h"
#include "measure/measure_handler.h"
#include "model/model_base.h"
#include "utils/numerical_stable.hpp"
#include "random.h"


namespace QuantumMonteCarlo {

    // alias conventions
    using TimeIndex = int;
    using RealScalar = double;
    using RealScalarVec = std::vector<double>;
    using ptrRealScalarVec = std::unique_ptr<std::vector<double>>;
    using SvdStack = Utils::SvdStack;
    using ptrSvdStack = std::unique_ptr<SvdStack>;
    using Matrix = Eigen::MatrixXd;
    using NumericalStable = Utils::NumericalStable;
    using GreensFunc = Eigen::MatrixXd;
    using GreensFuncVec = std::vector<Eigen::MatrixXd>;


    void DqmcWalker::set_physical_params( RealScalar beta, int time_size ) 
    {
        assert( beta > 0.0 );
        this->m_beta = beta;
        this->m_time_size = time_size;
        this->m_time_interval = beta / time_size;
    }


    void DqmcWalker::set_stabilization_pace( int stabilization_pace ) 
    {
        assert( stabilization_pace > 0 );
        this->m_stabilization_pace = stabilization_pace;
    }


    void DqmcWalker::initial( const LatticeBase& lattice, const MeasureHandler& meas_handler ) 
    {
        this->m_space_size = lattice.TotalSiteNum();
        this->m_current_time_slice = 0;
        this->m_wrap_error = 0.0;
        
        this->m_is_equaltime = meas_handler.isEqualTime();
        this->m_is_dynamic = meas_handler.isDynamic();
    }


    void DqmcWalker::initial_svd_stacks( const LatticeBase& lattice, const ModelBase& model ) 
    {
        // initialize udv stacks for sweep use
        // sweep process will start from 0 to beta, so we initialize svd_stack_right here.
        // stabilize the process every stabilization_pace steps

        // allocate memory for SvdStack class
        this->m_svd_stack_left_up = std::make_unique<SvdStack>(this->m_space_size, this->m_time_size);
        this->m_svd_stack_left_dn = std::make_unique<SvdStack>(this->m_space_size, this->m_time_size);
        this->m_svd_stack_right_up = std::make_unique<SvdStack>(this->m_space_size, this->m_time_size);
        this->m_svd_stack_right_dn = std::make_unique<SvdStack>(this->m_space_size, this->m_time_size);

        Matrix tmp_stack_up = Matrix::Identity(this->m_space_size, this->m_space_size);
        Matrix tmp_stack_dn = Matrix::Identity(this->m_space_size, this->m_space_size);

        // initial svd stacks for sweeping usages
        for (auto t = this->m_time_size; t >= 1; --t) {
            model.mult_transB_from_left(tmp_stack_up, t, +1.0);
            model.mult_transB_from_left(tmp_stack_dn, t, -1.0);

            // stabilize every nwrap steps with svd decomposition
            if ( (t - 1) % this->m_stabilization_pace == 0 ) {
                this->m_svd_stack_right_up->push(tmp_stack_up);
                this->m_svd_stack_right_dn->push(tmp_stack_dn);
                tmp_stack_up = Matrix::Identity(this->m_space_size, this->m_space_size);
                tmp_stack_dn = Matrix::Identity(this->m_space_size, this->m_space_size);
            }
        }
    }


    void DqmcWalker::initial_greens_function() 
    {
        // first allocate memory for greens functions
        this->m_green_tt_up = std::make_unique<GreensFunc>(this->m_space_size, this->m_space_size);
        this->m_green_tt_dn = std::make_unique<GreensFunc>(this->m_space_size, this->m_space_size);

        if ( this->m_is_equaltime || this->m_is_dynamic ) {
            this->m_vec_green_tt_up = std::make_unique<GreensFuncVec>(this->m_time_size, GreensFunc(this->m_space_size, this->m_space_size));
            this->m_vec_green_tt_dn = std::make_unique<GreensFuncVec>(this->m_time_size, GreensFunc(this->m_space_size, this->m_space_size));
        }

        if ( this->m_is_dynamic ) {
            this->m_green_t0_up = std::make_unique<GreensFunc>(this->m_space_size, this->m_space_size);
            this->m_green_t0_dn = std::make_unique<GreensFunc>(this->m_space_size, this->m_space_size);
            this->m_green_0t_up = std::make_unique<GreensFunc>(this->m_space_size, this->m_space_size);
            this->m_green_0t_dn = std::make_unique<GreensFunc>(this->m_space_size, this->m_space_size);

            this->m_vec_green_t0_up = std::make_unique<GreensFuncVec>(this->m_time_size, GreensFunc(this->m_space_size, this->m_space_size));
            this->m_vec_green_t0_dn = std::make_unique<GreensFuncVec>(this->m_time_size, GreensFunc(this->m_space_size, this->m_space_size));
            this->m_vec_green_0t_up = std::make_unique<GreensFuncVec>(this->m_time_size, GreensFunc(this->m_space_size, this->m_space_size));
            this->m_vec_green_0t_dn = std::make_unique<GreensFuncVec>(this->m_time_size, GreensFunc(this->m_space_size, this->m_space_size));
        }

        // compute greens function at time slice t = 0
        // which corresponds to imaginary-time tau = beta
        NumericalStable::compute_equaltime_greens( *this->m_svd_stack_left_up, *this->m_svd_stack_right_up, *this->m_green_tt_up );
        NumericalStable::compute_equaltime_greens( *this->m_svd_stack_left_dn, *this->m_svd_stack_right_dn, *this->m_green_tt_dn );
    }


    void DqmcWalker::initial_config_sign() 
    {   
        // allocate memory for config sign vector
        // if equal-time measurements are to be performed
        if ( this->m_is_equaltime ) {
            this->m_vec_config_sign = std::make_unique<RealScalarVec>(this->m_time_size, 0.0);
        }

        // initialize the sign of the initial bosonic configurations
        this->m_config_sign = ( this->m_green_tt_up->determinant() * this->m_green_tt_dn->determinant() >= 0 )? +1.0 : -1.0;
    }



    /*
     *  Update the aux bosonic fields at space-time position (t,i) 
     *  for all i with Metropolis probability, and, if the update is accepted, 
     *  perform a in-place update of the green's functions.
     *  Record the updated green's function at the life-end of this function.
     */
    void DqmcWalker::metropolis_update( ModelBase& model, TimeIndex t )
    {
        assert( this->m_current_time_slice == t );
        assert( t >= 0 && t <= this->m_time_size );

        const int eff_t = (t == 0)? this->m_time_size-1 : t-1;
        for (auto i = 0; i < this->m_space_size; ++i) {
            
            // obtain the ratio of flipping the bosonic field at (i,l)
            const auto update_ratio = model.get_update_ratio( *this, eff_t, i );

            if ( std::bernoulli_distribution(std::min(1.0, std::abs(update_ratio)))(Utils::Random::Engine) )
            {   
                // if accepted
                // update the greens functions
                model.update_greens_function( *this, eff_t, i );

                // update the bosonic fields
                model.update_bosonic_field( eff_t, i );

                // keep track of sign problem
                this->m_config_sign = ( update_ratio >= 0 )? +this->m_config_sign : -this->m_config_sign;
            }
        }
    }



    /*
     *  Propagate the greens functions from the current time slice t
     *  forwards to the time slice t+1 according to
     *      G(t+1) = B(t+1) * G(t) * B(t+1)^-1
     *  for both spin-1/2 states. 
     *  The greens functions are changed in place.
     */
    void DqmcWalker::wrap_from_0_to_beta( const ModelBase& model, TimeIndex t )
    {
        assert( t >= 0 && t <= this->m_time_size );

        const int eff_t = ( t == this->m_time_size )? 1 : t+1;
        model.mult_B_from_left     ( *this->m_green_tt_up, eff_t, +1 );
        model.mult_invB_from_right ( *this->m_green_tt_up, eff_t, +1 );
        model.mult_B_from_left     ( *this->m_green_tt_dn, eff_t, -1 );
        model.mult_invB_from_right ( *this->m_green_tt_dn, eff_t, -1 );
    }



    /*
     *  Propagate the greens function from the current time slice t
     *  downwards to the time slice t-1 according to
     *      G(t-1) = B(t)^-1 * G(t) * B(t)
     *  for both spin-1/2 states. 
     *  The greens functions are changed in place.
     */
    void DqmcWalker::wrap_from_beta_to_0( const ModelBase& model, TimeIndex t )
    {
        assert( t >= 0 && t <= this->m_time_size );

        const int eff_t = ( t == 0 )? this->m_time_size : t;
        model.mult_B_from_right   ( *this->m_green_tt_up, eff_t, +1 );
        model.mult_invB_from_left ( *this->m_green_tt_up, eff_t, +1 );
        model.mult_B_from_right   ( *this->m_green_tt_dn, eff_t, -1 );
        model.mult_invB_from_left ( *this->m_green_tt_dn, eff_t, -1 );
    }



    /*
     *  Update the space-time lattice of the auxiliary bosonic fields.
     *  For t = 1,2...,ts , attempt to update fields and propagate the greens functions
     *  Perform the stabilization every 'stabilization_pace' time slices
     */
    void DqmcWalker::sweep_from_0_to_beta( ModelBase& model )
    {
        this->m_current_time_slice++;

        const int stack_length = ( this->m_time_size % this->m_stabilization_pace == 0 )? 
                                   this->m_time_size/this->m_stabilization_pace 
                                 : this->m_time_size/this->m_stabilization_pace + 1 ;
        assert( this->m_current_time_slice == 1 );
        assert( this->m_svd_stack_left_up->empty() && this->m_svd_stack_left_dn->empty() );
        assert( this->m_svd_stack_right_up->StackLength() == stack_length && 
                this->m_svd_stack_right_dn->StackLength() == stack_length );

        // temporary matrices
        Matrix tmp_mat_up = Matrix::Identity(this->m_space_size, this->m_space_size);
        Matrix tmp_mat_dn = Matrix::Identity(this->m_space_size, this->m_space_size);

        // sweep upwards from 0 to beta
        for (auto t = 1; t <= this->m_time_size; ++t) 
        {
            // wrap green function to current time slice t
            this->wrap_from_0_to_beta( model, t-1 );

            // update auxiliary fields and record the updated greens functions
            this->metropolis_update( model, t );
            if ( this->m_is_equaltime ) {
                (*this->m_vec_green_tt_up)[t-1] = *this->m_green_tt_up;
                (*this->m_vec_green_tt_dn)[t-1] = *this->m_green_tt_dn;
                (*this->m_vec_config_sign)[t-1] = this->m_config_sign;
            }

            model.mult_B_from_left(tmp_mat_up, t, +1);
            model.mult_B_from_left(tmp_mat_dn, t, -1);

            // perform the stabilizations
            if ( t % this->m_stabilization_pace == 0 || t == this->m_time_size ) {
                // update svd stacks
                this->m_svd_stack_right_up->pop();
                this->m_svd_stack_right_dn->pop();
                this->m_svd_stack_left_up->push(tmp_mat_up);
                this->m_svd_stack_left_dn->push(tmp_mat_dn);

                // collect the wrapping errors
                Matrix tmp_green_tt_up = Matrix::Zero(this->m_space_size, this->m_space_size);
                Matrix tmp_green_tt_dn = Matrix::Zero(this->m_space_size, this->m_space_size);
                double tmp_wrap_error_tt_up = 0.0;
                double tmp_wrap_error_tt_dn = 0.0;

                // compute fresh greens every 'stabilization_pace' steps: g = ( 1 + stack_left*stack_right^T )^-1
                // stack_left = B(t-1) * ... * B(0)
                // stack_right = B(t)^T * ... * B(ts-1)^T
                NumericalStable::compute_equaltime_greens(*this->m_svd_stack_left_up, *this->m_svd_stack_right_up, tmp_green_tt_up);
                NumericalStable::compute_equaltime_greens(*this->m_svd_stack_left_dn, *this->m_svd_stack_right_dn, tmp_green_tt_dn);

                // compute wrapping errors
                NumericalStable::matrix_compare_error(tmp_green_tt_up, *this->m_green_tt_up, tmp_wrap_error_tt_up);
                NumericalStable::matrix_compare_error(tmp_green_tt_dn, *this->m_green_tt_dn, tmp_wrap_error_tt_dn);
                this->m_wrap_error = std::max(this->m_wrap_error, std::max(tmp_wrap_error_tt_up, tmp_wrap_error_tt_dn));

                *this->m_green_tt_up = tmp_green_tt_up;
                *this->m_green_tt_dn = tmp_green_tt_dn;

                if ( this->m_is_equaltime ) {
                    (*this->m_vec_green_tt_up)[t-1] = *this->m_green_tt_up;
                    (*this->m_vec_green_tt_dn)[t-1] = *this->m_green_tt_dn;
                }

                tmp_mat_up = Matrix::Identity(this->m_space_size, this->m_space_size);
                tmp_mat_dn = Matrix::Identity(this->m_space_size, this->m_space_size);
            }

            // finally stop at time slice t = ts + 1
            this->m_current_time_slice++;
        }

        // end with fresh greens functions
        if ( this->m_is_equaltime ) {
            (*this->m_vec_green_tt_up)[this->m_time_size-1] = *this->m_green_tt_up;
            (*this->m_vec_green_tt_dn)[this->m_time_size-1] = *this->m_green_tt_dn;
        }
    }
    


    /*
     *  Update the space-time lattice of the auxiliary bosonic fields.
     *  For l = ts,ts-1,...,1 , attempt to update fields and propagate the greens functions
     *  Perform the stabilization every 'stabilization_pace' time slices
     */
    void DqmcWalker::sweep_from_beta_to_0( ModelBase& model )
    {
        this->m_current_time_slice--;

        const int stack_length = ( this->m_time_size % this->m_stabilization_pace == 0 )? 
                                   this->m_time_size/this->m_stabilization_pace 
                                 : this->m_time_size/this->m_stabilization_pace + 1;
        assert( this->m_current_time_slice == this->m_time_size );
        assert( this->m_svd_stack_right_up->empty() && this->m_svd_stack_right_dn->empty() );
        assert( this->m_svd_stack_left_up->StackLength() == stack_length && 
                this->m_svd_stack_left_dn->StackLength() == stack_length );

        // temporary matrices
        Matrix tmp_mat_up = Matrix::Identity(this->m_space_size, this->m_space_size);
        Matrix tmp_mat_dn = Matrix::Identity(this->m_space_size, this->m_space_size);

        // sweep downwards from beta to 0
        for (auto t = this->m_time_size; t >= 1; --t) {

            // perform the stabilizations
            if ( t % this->m_stabilization_pace == 0 && t != this->m_time_size ) {
                // update svd stacks
                this->m_svd_stack_left_up->pop();
                this->m_svd_stack_left_dn->pop();
                this->m_svd_stack_right_up->push(tmp_mat_up);
                this->m_svd_stack_right_dn->push(tmp_mat_dn);

                // collect the wrapping errors
                Matrix tmp_green_tt_up = Matrix::Zero(this->m_space_size, this->m_space_size);
                Matrix tmp_green_tt_dn = Matrix::Zero(this->m_space_size, this->m_space_size);
                double tmp_wrap_error_tt_up = 0.0;
                double tmp_wrap_error_tt_dn = 0.0;

                NumericalStable::compute_equaltime_greens(*this->m_svd_stack_left_up, *this->m_svd_stack_right_up, tmp_green_tt_up);
                NumericalStable::compute_equaltime_greens(*this->m_svd_stack_left_dn, *this->m_svd_stack_right_dn, tmp_green_tt_dn);

                // compute the wrapping errors
                NumericalStable::matrix_compare_error(tmp_green_tt_up, *this->m_green_tt_up, tmp_wrap_error_tt_up);
                NumericalStable::matrix_compare_error(tmp_green_tt_dn, *this->m_green_tt_dn, tmp_wrap_error_tt_dn);
                this->m_wrap_error = std::max(this->m_wrap_error, std::max(tmp_wrap_error_tt_up, tmp_wrap_error_tt_dn));

                *this->m_green_tt_up = tmp_green_tt_up;
                *this->m_green_tt_dn = tmp_green_tt_dn;

                tmp_mat_up = Matrix::Identity(this->m_space_size, this->m_space_size);
                tmp_mat_dn = Matrix::Identity(this->m_space_size, this->m_space_size);
            }

            // update auxiliary fields and record the updated greens functions
            this->metropolis_update( model, t );
            if ( this->m_is_equaltime ) {
                (*this->m_vec_green_tt_up)[t-1] = *this->m_green_tt_up;
                (*this->m_vec_green_tt_dn)[t-1] = *this->m_green_tt_dn;
                (*this->m_vec_config_sign)[t-1] = this->m_config_sign;
            }

            model.mult_transB_from_left(tmp_mat_up, t, +1);
            model.mult_transB_from_left(tmp_mat_dn, t, -1);

            this->wrap_from_beta_to_0( model, t );

            this->m_current_time_slice--;
        }

        // at time slice t = 0
        this->m_svd_stack_left_up->pop();
        this->m_svd_stack_left_dn->pop();
        this->m_svd_stack_right_up->push(tmp_mat_up);
        this->m_svd_stack_right_dn->push(tmp_mat_dn);

        NumericalStable::compute_equaltime_greens(*this->m_svd_stack_left_up, *this->m_svd_stack_right_up, *this->m_green_tt_up);
        NumericalStable::compute_equaltime_greens(*this->m_svd_stack_left_dn, *this->m_svd_stack_right_dn, *this->m_green_tt_dn);

        // end with fresh greens functions
        if ( this->m_is_equaltime ) {
            (*this->m_vec_green_tt_up)[this->m_time_size-1] = *this->m_green_tt_up;
            (*this->m_vec_green_tt_dn)[this->m_time_size-1] = *this->m_green_tt_dn;
        }
    }



    /*
     *  Calculate time-displaced (dynamical) greens functions, while the auxiliary fields remain unchanged.
     *  For l = 1,2...,ts , recompute the SvdStacks every 'stabilization_pace' time slices.
     *  The collected dynamic greens functions are stored in m_vec_green_t0(0t)_up(dn).
     *  Note that the equal-time greens functions are also re-calculated 
     *  according to the current auxiliary field configurations, 
     *  which are stored in m_vec_green_tt_up(dn).
     */
    void DqmcWalker::sweep_for_dynamic_greens( ModelBase& model )
    {
        if ( this->m_is_dynamic ) {

            this->m_current_time_slice++;
            const int stack_length = ( this->m_time_size % this->m_stabilization_pace == 0 )? 
                                       this->m_time_size/this->m_stabilization_pace 
                                     : this->m_time_size/this->m_stabilization_pace + 1 ;
            assert( this->m_current_time_slice == 1 );
            assert( this->m_svd_stack_left_up->empty() && this->m_svd_stack_left_dn->empty() );
            assert( this->m_svd_stack_right_up->StackLength() == stack_length && 
                    this->m_svd_stack_right_dn->StackLength() == stack_length );

            // initialize greens functions: at t = 0, gt0 = g00, g0t = g00 - 1
            *this->m_green_t0_up = *this->m_green_tt_up;
            *this->m_green_t0_dn = *this->m_green_tt_dn;
            *this->m_green_0t_up = *this->m_green_tt_up - Matrix::Identity(this->m_space_size, this->m_space_size);
            *this->m_green_0t_dn = *this->m_green_tt_dn - Matrix::Identity(this->m_space_size, this->m_space_size);

            // temporary matrices
            Matrix tmp_mat_up = Matrix::Identity(this->m_space_size, this->m_space_size);
            Matrix tmp_mat_dn = Matrix::Identity(this->m_space_size, this->m_space_size);

            // sweep forwards from 0 to beta
            for (auto t = 1; t <= this->m_time_size; ++t) {
                // wrap the equal time greens functions to current time slice t
                this->wrap_from_0_to_beta( model, t-1 );
                (*this->m_vec_green_tt_up)[t-1] = *this->m_green_tt_up;
                (*this->m_vec_green_tt_dn)[t-1] = *this->m_green_tt_dn;
            
                // calculate and record the time-displaced greens functions at different time slices
                model.mult_B_from_left(*this->m_green_t0_up, t, +1);
                model.mult_B_from_left(*this->m_green_t0_dn, t, -1);
                (*this->m_vec_green_t0_up)[t-1] = *this->m_green_t0_up;
                (*this->m_vec_green_t0_dn)[t-1] = *this->m_green_t0_dn;

                model.mult_invB_from_right(*this->m_green_0t_up, t, +1);
                model.mult_invB_from_right(*this->m_green_0t_dn, t, -1);
                (*this->m_vec_green_0t_up)[t-1] = *this->m_green_0t_up;
                (*this->m_vec_green_0t_dn)[t-1] = *this->m_green_0t_dn;

                model.mult_B_from_left(tmp_mat_up, t, +1);
                model.mult_B_from_left(tmp_mat_dn, t, -1);

                // perform the stabilizations
                if ( t % this->m_stabilization_pace == 0 || t == this->m_time_size ) {
                    // update svd stacks
                    this->m_svd_stack_right_up->pop();
                    this->m_svd_stack_right_dn->pop();
                    this->m_svd_stack_left_up->push(tmp_mat_up);
                    this->m_svd_stack_left_dn->push(tmp_mat_dn);

                    // collect the wrapping errors
                    Matrix tmp_green_t0_up = Matrix::Zero(this->m_space_size, this->m_space_size);
                    Matrix tmp_green_t0_dn = Matrix::Zero(this->m_space_size, this->m_space_size);
                    Matrix tmp_green_0t_up = Matrix::Zero(this->m_space_size, this->m_space_size);
                    Matrix tmp_green_0t_dn = Matrix::Zero(this->m_space_size, this->m_space_size);
                    double tmp_wrap_error_t0_up = 0.0;
                    double tmp_wrap_error_t0_dn = 0.0;
                    double tmp_wrap_error_0t_up = 0.0;
                    double tmp_wrap_error_0t_dn = 0.0;

                    // compute fresh greens every nwrap steps
                    // stack_left = B(t-1) * ... * B(0)
                    // stack_right = B(t)^T * ... * B(ts-1)^T
                    // equal time green's function are re-evaluated for current field configurations
                    NumericalStable::compute_equaltime_greens (*this->m_svd_stack_left_up, *this->m_svd_stack_right_up, *this->m_green_tt_up);
                    NumericalStable::compute_equaltime_greens (*this->m_svd_stack_left_dn, *this->m_svd_stack_right_dn, *this->m_green_tt_dn);
                    NumericalStable::compute_dynamic_greens   (*this->m_svd_stack_left_up, *this->m_svd_stack_right_up, tmp_green_t0_up, tmp_green_0t_up);
                    NumericalStable::compute_dynamic_greens   (*this->m_svd_stack_left_dn, *this->m_svd_stack_right_dn, tmp_green_t0_dn, tmp_green_0t_dn);

                    // compute wrapping errors
                    NumericalStable::matrix_compare_error(tmp_green_t0_up, *this->m_green_t0_up, tmp_wrap_error_t0_up);
                    NumericalStable::matrix_compare_error(tmp_green_t0_dn, *this->m_green_t0_dn, tmp_wrap_error_t0_dn);
                    this->m_wrap_error = std::max(this->m_wrap_error, std::max(tmp_wrap_error_t0_up, tmp_wrap_error_t0_dn));

                    NumericalStable::matrix_compare_error(tmp_green_0t_up, *this->m_green_0t_up, tmp_wrap_error_0t_up);
                    NumericalStable::matrix_compare_error(tmp_green_0t_dn, *this->m_green_0t_dn, tmp_wrap_error_0t_dn);
                    this->m_wrap_error = std::max(this->m_wrap_error, std::max(tmp_wrap_error_0t_up, tmp_wrap_error_0t_dn));

                    *this->m_green_t0_up = tmp_green_t0_up;
                    *this->m_green_t0_dn = tmp_green_t0_dn;
                    *this->m_green_0t_up = tmp_green_0t_up;
                    *this->m_green_0t_dn = tmp_green_0t_dn;

                    (*this->m_vec_green_tt_up)[t-1] = *this->m_green_tt_up;
                    (*this->m_vec_green_tt_dn)[t-1] = *this->m_green_tt_dn;
                    (*this->m_vec_green_t0_up)[t-1] = *this->m_green_t0_up;
                    (*this->m_vec_green_t0_dn)[t-1] = *this->m_green_t0_dn;
                    (*this->m_vec_green_0t_up)[t-1] = *this->m_green_0t_up;
                    (*this->m_vec_green_0t_dn)[t-1] = *this->m_green_0t_dn;

                    tmp_mat_up = Matrix::Identity(this->m_space_size, this->m_space_size);
                    tmp_mat_dn = Matrix::Identity(this->m_space_size, this->m_space_size);
                }

                // finally stop at time slice t = ts + 1
                this->m_current_time_slice++;
            }
        }
    }


} // namespace QuantumMonteCarlo