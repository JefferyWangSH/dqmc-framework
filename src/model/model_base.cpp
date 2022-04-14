#include "model/model_base.h"


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



} // namespace Model