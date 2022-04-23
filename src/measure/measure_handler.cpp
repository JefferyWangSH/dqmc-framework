#include "measure/measure_handler.h"
#include "model/model_base.h"
#include "lattice/lattice_base.h"


namespace Measure {

    const bool MeasureHandler::isWarmUp() const { return this->m_is_warmup; }
    const bool MeasureHandler::isEqualTime() const { return this->m_is_equaltime; }
    const bool MeasureHandler::isDynamic() const { return this->m_is_dynamic; }


    void MeasureHandler::set_measure_params( int sweeps_warmup, int bin_num, int bin_size, int sweeps_between_bins )
    {
        assert( sweeps_warmup >= 0 );
        assert( bin_num >= 0 );
        assert( bin_size >= 0 );
        assert( sweeps_warmup >= 0 );
        this->m_sweeps_warmup = sweeps_warmup;
        this->m_bin_num = bin_num;
        this->m_bin_size = bin_size;
        this->m_sweeps_between_bins = sweeps_between_bins;
    }


    void MeasureHandler::set_observables( ObsList obs_list )
    {
        this->m_obs_list = obs_list;
    }


    void MeasureHandler::initial()
    {
        // initialize ObservableHandler
        Observable::ObservableHandler::initial(this->m_obs_list);

        this->m_is_warmup = (this->m_sweeps_warmup != 0);
        this->m_is_equaltime = (   !this->m_eqtime_scalar_obs.empty() 
                                || !this->m_eqtime_vector_obs.empty() 
                                || !this->m_eqtime_matrix_obs.empty()  );
        this->m_is_dynamic   = (   !this->m_dynamic_scalar_obs.empty() 
                                || !this->m_dynamic_vector_obs.empty() 
                                || !this->m_dynamic_matrix_obs.empty() );
    }


    


} // namespace Measure
