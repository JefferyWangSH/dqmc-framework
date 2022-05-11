#include "lattice/lattice_base.h"

namespace Lattice {
    
    const LatticeBool LatticeBase::InitialStatus() const { return this->m_initial_status; }
    const LatticeInt  LatticeBase::SpaceDim() const { return this->m_space_dim; }
    const LatticeInt  LatticeBase::SpaceSize() const { return this->m_space_size; }
    const LatticeInt  LatticeBase::SideLength() const { return this->m_side_length; }
    const LatticeInt  LatticeBase::CoordinationNumber() const { return this->m_coordination_number; }
    const LatticeInt  LatticeBase::kStarsNum() const { return this->m_num_k_stars; }


    const LatticeInt LatticeBase::GammaPointIndex() const { return this->m_gamma_point_index; }
    const LatticeInt LatticeBase::MPointIndex() const { return this->m_m_point_index; }
    const LatticeInt LatticeBase::XPointIndex() const { return this->m_x_point_index; }
    const LatticeIntVec& LatticeBase::DeltaLineIndex() const { return this->m_delta_line_index; }
    const LatticeIntVec& LatticeBase::ZLineIndex() const { return this->m_z_line_index; }
    const LatticeIntVec& LatticeBase::SigmaLineIndex() const { return this->m_sigma_line_index; }
    const LatticeIntVec& LatticeBase::kStarsIndex() const { return this->m_k_stars_index; }
    const LatticeIntVec& LatticeBase::Gamma2X2M2GammaLoopIndex() const { return this->m_gamma2x2m2gamma_loop_index; }


    const MatrixDouble& LatticeBase::HoppingMatrix() const { return this->m_hopping_matrix; }
    const MatrixDouble& LatticeBase::FourierFactor() const { return this->m_fourier_factor_table; }


    const LatticeInt LatticeBase::NearestNeighbour( const LatticeInt site_index, const LatticeInt direction ) const
    {
        assert( site_index >= 0 && site_index < this->m_space_size );
        assert( direction >= 0 && direction < this->m_coordination_number );
        return this->m_nearest_neighbour_table(site_index, direction);
    }


    const LatticeDouble LatticeBase::FourierFactor( const LatticeInt site_index, const LatticeInt momentum_index ) const
    {
        assert( site_index >= 0 && site_index < this->m_space_size );
        assert( momentum_index >= 0 && momentum_index < this->m_num_k_stars );
        return this->m_fourier_factor_table(site_index, momentum_index);
    }


    const VectorInt LatticeBase::Index2Site( const LatticeInt site_index ) const 
    {   
        assert( site_index >= 0 && site_index < this->m_space_size );
        return this->m_index2site_table.row(site_index);
    }


    const LatticeInt LatticeBase::Index2Site( const LatticeInt site_index, const LatticeInt axis ) const 
    {   
        assert( site_index >= 0 && site_index < this->m_space_size );
        assert( axis >= 0 && axis < this->m_space_dim );
        return this->m_index2site_table(site_index, axis);
    }


    const VectorDouble LatticeBase::Index2Momentum( const LatticeInt momentum_index ) const 
    {
        assert( momentum_index >= 0 && momentum_index < this->m_num_k_stars );
        return this->m_index2momentum_table.row(momentum_index);
    }


    const LatticeDouble LatticeBase::Index2Momentum( const LatticeInt momentum_index, const LatticeInt axis ) const 
    {   
        assert( momentum_index >= 0 && momentum_index < this->m_num_k_stars );
        assert( axis >= 0 && axis < this->m_space_dim );
        return this->m_index2momentum_table(momentum_index, axis);
    }


    const LatticeInt LatticeBase::Displacement( const LatticeInt site1_index, const LatticeInt site2_index ) const
    {
        assert( site1_index >= 0 && site1_index < this->m_space_size );
        assert( site2_index >= 0 && site2_index < this->m_space_size );
        return this->m_displacement_table(site1_index, site2_index);
    }


} // namespace Lattice
