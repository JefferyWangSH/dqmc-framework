#include "dqmc_walker.h"
#include "svd_stack.h"
#include "lattice/lattice_base.h"
#include "measure/measure_handler.h"
#include "model/model_base.h"
#include "utils/numerical_stable.hpp"


namespace QuantumMonteCarlo {

    // alias conventions
    using TimeIndex = int;
    using RealScalar = double;
    using ptrRealScalarVec = std::unique_ptr<std::vector<double>>;
    using SvdStack = Utils::SvdStack;
    using ptrSvdStack = std::unique_ptr<SvdStack>;
    using Matrix = Eigen::MatrixXd;
    using NumericalStable = Utils::NumericalStable;


    const int DqmcWalker::TimeSliceNum() const {
        return this->m_time_size;
    }

    const RealScalar DqmcWalker::Beta() const {
        return this->m_beta;
    }

    const RealScalar DqmcWalker::TimeInterval() const {
        return this->m_time_interval;
    }

    const RealScalar DqmcWalker::WrapError() const {
        return this->m_wrap_error;
    }

    const int DqmcWalker::StabilizationPace() const {
        return this->m_stabilization_pace;
    }


    void DqmcWalker::set_physical_params( RealScalar beta, int time_size ) {
        assert( beta > 0.0 );
        this->m_beta = beta;
        this->m_time_size = time_size;
        this->m_time_interval = beta / time_size;
    }


    void DqmcWalker::set_stabilization_pace( int stabilization_pace ) {
        assert( stabilization_pace > 0 );
        this->m_stabilization_pace = stabilization_pace;
    }


    void DqmcWalker::initial( const LatticeBase& lattice, const MeasureHandler& meas_handler ) {
        this->m_space_size = lattice.TotalSiteNum();
        this->m_current_time_slice = 0;
        this->m_wrap_error = 0.0;
        
        this->m_is_equaltime = meas_handler.isEqualTime();
        this->m_is_dynamic = meas_handler.isDynamic();
    }


    void DqmcWalker::initial_svd_stacks( const LatticeBase& lattice, const ModelBase& model ) {
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


    void DqmcWalker::initial_greens_function( ModelBase& model ) {
        // initialize greens function at time slice t = 0
        // which corresponds to imaginary-time tau = beta
        NumericalStable::compute_greens_eqtime( *this->m_svd_stack_left_up, *this->m_svd_stack_right_up, model.GreenttUp() );
        NumericalStable::compute_greens_eqtime( *this->m_svd_stack_left_dn, *this->m_svd_stack_right_dn, model.GreenttDn() );
    }




} // namespace QuantumMonteCarlo