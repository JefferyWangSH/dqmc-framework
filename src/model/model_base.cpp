#include "model/model_base.h"
#include "checkerboard/checkerboard_base.h"


namespace Model {

    void ModelBase::mult_expK_from_left( GreensFunc& green ) const 
    { 
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        green = this->m_expK_mat * green;
    }

    void ModelBase::mult_expK_from_right( GreensFunc& green ) const 
    { 
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        green = green * this->m_expK_mat; 
    }

    void ModelBase::mult_inv_expK_from_left( GreensFunc& green ) const 
    { 
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        green = this->m_inv_expK_mat * green; 
    }

    void ModelBase::mult_inv_expK_from_right( GreensFunc& green ) const 
    { 
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        green = green * this->m_inv_expK_mat; 
    }
    
    void ModelBase::mult_trans_expK_from_left( GreensFunc& green ) const 
    { 
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        green = this->m_trans_expK_mat * green;
    }


    void ModelBase::link()
    {
        this->m_mult_expK_from_left       = std::bind(&ModelBase::mult_expK_from_left, this, std::placeholders::_1);
        this->m_mult_expK_from_right      = std::bind(&ModelBase::mult_expK_from_right, this, std::placeholders::_1);
        this->m_mult_inv_expK_from_left   = std::bind(&ModelBase::mult_inv_expK_from_left, this, std::placeholders::_1);
        this->m_mult_inv_expK_from_right  = std::bind(&ModelBase::mult_inv_expK_from_right, this, std::placeholders::_1);
        this->m_mult_trans_expK_from_left = std::bind(&ModelBase::mult_trans_expK_from_left, this, std::placeholders::_1);
    }


    void ModelBase::link( const CheckerBoardBase& checkerboard )
    {
        this->m_mult_expK_from_left       = std::bind(&CheckerBoardBase::mult_expK_from_left, &checkerboard, std::placeholders::_1);
        this->m_mult_expK_from_right      = std::bind(&CheckerBoardBase::mult_expK_from_right, &checkerboard, std::placeholders::_1);
        this->m_mult_inv_expK_from_left   = std::bind(&CheckerBoardBase::mult_inv_expK_from_left, &checkerboard, std::placeholders::_1);
        this->m_mult_inv_expK_from_right  = std::bind(&CheckerBoardBase::mult_inv_expK_from_right, &checkerboard, std::placeholders::_1);
        this->m_mult_trans_expK_from_left = std::bind(&CheckerBoardBase::mult_trans_expK_from_left, &checkerboard, std::placeholders::_1);
    }


} // namespace Model