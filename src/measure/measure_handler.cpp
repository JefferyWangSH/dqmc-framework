#include "measure/measure_handler.h"
#include "model/model_base.h"
#include "lattice/lattice_base.h"
#include "dqmc_walker.h"


namespace Measure {

    const bool MeasureHandler::isWarmUp() const { return this->m_is_warmup; }
    const bool MeasureHandler::isEqualTime() const { return this->m_is_equaltime; }
    const bool MeasureHandler::isDynamic() const { return this->m_is_dynamic; }

    const int MeasureHandler::WarmUpSweeps() const { return this->m_sweeps_warmup; }
    const int MeasureHandler::SweepsBetweenBins() const { return this->m_sweeps_between_bins; }
    const int MeasureHandler::BinsNum() const { return this->m_bin_num; }
    const int MeasureHandler::BinsSize() const { return this->m_bin_size; }


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


    void MeasureHandler::initial( const LatticeBase& lattice, const DqmcWalker& walker )
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

        // set up parameters for the observables
        // for equal-time observables
        if ( this->m_is_equaltime ) {
            // allocate for equal-time sign measurements
            this->m_equaltime_sign->set_zero_element(0.0);
            this->m_equaltime_sign->set_size_of_bin(this->m_bin_size);
            this->m_equaltime_sign->allocate();

            for (auto& scalar_obs : this->m_eqtime_scalar_obs) {
                scalar_obs->set_zero_element(0.0);
                scalar_obs->set_size_of_bin(this->m_bin_size);
                scalar_obs->allocate();
            }
            for (auto& vector_obs : this->m_eqtime_vector_obs) {
                vector_obs->set_zero_element(Vector::Zero(walker.TimeSize()));
                vector_obs->set_size_of_bin(this->m_bin_size);
                vector_obs->allocate();
            }
            for (auto& matrix_obs : this->m_eqtime_matrix_obs) {
                // the dimensions of the observable can be adjusted or specialized
                // todo
                matrix_obs->set_zero_element(Matrix::Zero(lattice.SpaceSize(), lattice.SpaceSize()));
                matrix_obs->set_size_of_bin(this->m_bin_size);
                matrix_obs->allocate();
            }
        }

        // for dynamic observables
        if ( this->m_is_dynamic ) {
            // allocate for dynamic sign measurements
            this->m_dynamic_sign->set_zero_element(0.0);
            this->m_dynamic_sign->set_size_of_bin(this->m_bin_size);
            this->m_dynamic_sign->allocate();

            for (auto& scalar_obs : this->m_dynamic_scalar_obs) {
                scalar_obs->set_zero_element(0.0);
                scalar_obs->set_size_of_bin(this->m_bin_size);
                scalar_obs->allocate();
            }
            for (auto& vector_obs : this->m_dynamic_vector_obs) {
                vector_obs->set_zero_element(Vector::Zero(walker.TimeSize()));
                vector_obs->set_size_of_bin(this->m_bin_size);
                vector_obs->allocate();
            }
            for (auto& matrix_obs : this->m_dynamic_matrix_obs) {
                // the dimensions of the observable can be adjusted or specialized
                // todo
                matrix_obs->set_zero_element(Matrix::Zero(lattice.SpaceSize(), lattice.SpaceSize()));
                matrix_obs->set_size_of_bin(this->m_bin_size);
                matrix_obs->allocate();
            }
        }
    }


    void MeasureHandler::equaltime_measure( const ModelBase& model, const LatticeBase& lattice )
    {
        for (auto& scalar_obs : this->m_eqtime_scalar_obs) {
            scalar_obs->measure(*this, model, lattice);
        }
        for (auto& vector_obs : this->m_eqtime_vector_obs) {
            vector_obs->measure(*this, model, lattice);
        }
        for (auto& matrix_obs : this->m_eqtime_matrix_obs) {
            matrix_obs->measure(*this, model, lattice);
        }
        this->m_equaltime_sign->measure(*this, model, lattice);
    }

    
    void MeasureHandler::dynamic_measure( const ModelBase& model, const LatticeBase& lattice )
    {
        for (auto& scalar_obs : this->m_dynamic_scalar_obs) {
            scalar_obs->measure(*this, model, lattice);
        }
        for (auto& vector_obs : this->m_dynamic_vector_obs) {
            vector_obs->measure(*this, model, lattice);
        }
        for (auto& matrix_obs : this->m_dynamic_matrix_obs) {
            matrix_obs->measure(*this, model, lattice);
        }
        this->m_dynamic_sign->measure(*this, model, lattice);
    }


    void MeasureHandler::normalize_stats()
    {   
        if ( this->m_is_equaltime ) {
            // normalize the sign measurment first
            this->m_equaltime_sign->tmp_value() /= this->m_equaltime_sign->counts();

            // normalize observables by the countings and the mean value of the sign
            for (auto& scalar_obs : this->m_eqtime_scalar_obs) {
                scalar_obs->tmp_value() /= scalar_obs->counts() * this->m_equaltime_sign->tmp_value();
            }
            for (auto& vector_obs : this->m_eqtime_vector_obs) {
                vector_obs->tmp_value() /= vector_obs->counts() * this->m_equaltime_sign->tmp_value();
            }
            for (auto& matrix_obs : this->m_eqtime_matrix_obs) {
                matrix_obs->tmp_value() /= matrix_obs->counts() * this->m_equaltime_sign->tmp_value();
            }
        }

        if ( this->m_is_dynamic ) {
            // normalize the sign measurment first
            this->m_dynamic_sign->tmp_value() /= this->m_dynamic_sign->counts();

            for (auto& scalar_obs : this->m_dynamic_scalar_obs) {
                scalar_obs->tmp_value() /= scalar_obs->counts() * this->m_dynamic_sign->tmp_value();
            }
            for (auto& vector_obs : this->m_dynamic_vector_obs) {
                vector_obs->tmp_value() /= vector_obs->counts() * this->m_dynamic_sign->tmp_value();
            }
            for (auto& matrix_obs : this->m_dynamic_matrix_obs) {
                matrix_obs->tmp_value() /= matrix_obs->counts() * this->m_dynamic_sign->tmp_value();  
            }
        }
    }


    void MeasureHandler::write_stats_to_bins( int bin )
    {
        if ( this->m_is_equaltime ) {
            this->m_equaltime_sign->bin_data(bin) = this->m_equaltime_sign->tmp_value();

            for (auto& scalar_obs : this->m_eqtime_scalar_obs) {
                scalar_obs->bin_data(bin) = scalar_obs->tmp_value();
            }
            for (auto& vector_obs : this->m_eqtime_vector_obs) {
                vector_obs->bin_data(bin) = vector_obs->tmp_value();
            }
            for (auto& matrix_obs : this->m_eqtime_matrix_obs) {
                matrix_obs->bin_data(bin) = matrix_obs->tmp_value();
            }
        }

        if ( this->m_is_dynamic ) {
            this->m_dynamic_sign->bin_data(bin) = this->m_dynamic_sign->tmp_value();

            for (auto& scalar_obs : this->m_dynamic_scalar_obs) {
                scalar_obs->bin_data(bin) = scalar_obs->tmp_value();
            }
            for (auto& vector_obs : this->m_dynamic_vector_obs) {
                vector_obs->bin_data(bin) = vector_obs->tmp_value();
            }
            for (auto& matrix_obs : this->m_dynamic_matrix_obs) {
                matrix_obs->bin_data(bin) = matrix_obs->tmp_value();
            }
        }
    }


    void MeasureHandler::analyse_stats()
    {
        if ( this->m_is_equaltime ) {
            this->m_equaltime_sign->analyse();
            for (auto& scalar_obs : this->m_eqtime_scalar_obs) { scalar_obs->analyse(); }
            for (auto& vector_obs : this->m_eqtime_vector_obs) { vector_obs->analyse(); }
            for (auto& matrix_obs : this->m_eqtime_matrix_obs) { matrix_obs->analyse(); }
        }

        if ( this->m_is_dynamic ) {
            this->m_dynamic_sign->analyse();
            for (auto& scalar_obs : this->m_dynamic_scalar_obs) { scalar_obs->analyse(); }
            for (auto& vector_obs : this->m_dynamic_vector_obs) { vector_obs->analyse(); }
            for (auto& matrix_obs : this->m_dynamic_matrix_obs) { matrix_obs->analyse(); }
        }
    }
    
    
    void MeasureHandler::clear_temporary()
    {   
        // clear the temporary data for all the observables
        if ( this->m_is_equaltime ) {
            this->m_equaltime_sign->clear_temporary();
            for (auto& scalar_obs : this->m_eqtime_scalar_obs) { scalar_obs->clear_temporary(); }
            for (auto& vector_obs : this->m_eqtime_vector_obs) { vector_obs->clear_temporary(); }
            for (auto& matrix_obs : this->m_eqtime_matrix_obs) { matrix_obs->clear_temporary(); }
        }

        if ( this->m_is_dynamic ) {
            this->m_dynamic_sign->clear_temporary();
            for (auto& scalar_obs : this->m_dynamic_scalar_obs) { scalar_obs->clear_temporary(); }
            for (auto& vector_obs : this->m_dynamic_vector_obs) { vector_obs->clear_temporary(); }
            for (auto& matrix_obs : this->m_dynamic_matrix_obs) { matrix_obs->clear_temporary(); }
        }
    }


} // namespace Measure
