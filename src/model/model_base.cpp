#include "model/model_base.h"
#include "lattice/lattice_base.h"
#include "measure/measure_handler.h"
#include "dqmc_walker.h"


namespace Model {

    GreensFunc& ModelBase::GreenttUp() {
        // todo: this may cause problems if the pointer is nullptr
        return *this->m_green_tt_up;
    }

    GreensFunc& ModelBase::GreenttDn() {
        return *this->m_green_tt_dn;
    }

    GreensFunc& ModelBase::Greent0Up() {
        return *this->m_green_t0_up;
    }

    GreensFunc& ModelBase::Greent0Dn() {
        return *this->m_green_t0_dn;
    }

    GreensFunc& ModelBase::Green0tUp() {
        return *this->m_green_0t_up;
    }

    GreensFunc& ModelBase::Green0tDn() {
        return *this->m_green_0t_dn;
    }

    GreensFuncVec& ModelBase::vecGreenttUp() {
        return *this->m_vec_green_tt_up;
    }

    GreensFuncVec& ModelBase::vecGreenttDn() {
        return *this->m_vec_green_tt_dn;
    }

    GreensFuncVec& ModelBase::vecGreent0Up() {
        return *this->m_vec_green_t0_up;
    }

    GreensFuncVec& ModelBase::vecGreent0Dn() {
        return *this->m_vec_green_t0_dn;
    }

    GreensFuncVec& ModelBase::vecGreen0tUp() {
        return *this->m_vec_green_0t_up;
    }

    GreensFuncVec& ModelBase::vecGreen0tDn() {
        return *this->m_vec_green_0t_dn;
    }


    void ModelBase::initial_greens_function(    const Lattice& lattice, 
                                                const Walker& walker, 
                                                const MeasureHandler& meas_handler ) 
    {
        this->m_space_size = lattice.TotalSiteNum();
        this->m_time_size  = walker.TimeSliceNum();
        this->m_is_equaltime = meas_handler.isEqualTime();
        this->m_is_dynamic = meas_handler.isDynamic();
        
        // allocate memory for greens functions
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

    }



} // namespace Model