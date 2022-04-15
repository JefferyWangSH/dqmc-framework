#include "dqmc_walker.h"
#include "svd_stack.h"


namespace QuantumMonteCarlo {

    // alias conventions
    using TimeIndex = int;
    using RealScalar = double;
    using ptrRealScalarVec = std::unique_ptr<std::vector<double>>;
    using SvdStack = Utils::SvdStack;
    using ptrSvdStack = std::unique_ptr<SvdStack>;


    const int DqmcWalker::TimeSliceNum() const {
        return this->m_time_slice_num;
    }

    const RealScalar DqmcWalker::Beta() const {
        return this->m_beta;
    }

    const RealScalar DqmcWalker::TimeInterval() const {
        return this->m_time_interval;
    }


    void DqmcWalker::set_physical_params( RealScalar beta, int time_slice_num ) {
        assert( beta > 0.0 );
        this->m_beta = beta;
        this->m_time_slice_num = time_slice_num;
        this->m_time_interval = beta / time_slice_num;
    }



} // namespace QuantumMonteCarlo