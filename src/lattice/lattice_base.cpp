#include "lattice/lattice_base.h"

namespace Lattice {
    
    const LatticeInt LatticeBase::SpaceDim() const { return this->m_space_dim; }
    const LatticeInt LatticeBase::SpaceSize() const { return this->m_space_size; }
    const LatticeInt LatticeBase::SideLength() const { return this->m_side_length; }
    const MatrixDouble& LatticeBase::HoppingMatrix() const { return this->m_hopping_matrix; }


    const LatticeInt LatticeBase::NearestNeighbour(const LatticeInt site_index, const LatticeInt direction) const {
        return this->m_nearest_neighbour_table(site_index, direction);
    }

    const VectorInt LatticeBase::NearestNeighbour(const LatticeInt site_index) const {
        return this->m_nearest_neighbour_table.row(site_index);
    }

    const VectorInt LatticeBase::index2site(const LatticeInt site_index) const {
        return this->m_index2site_table.row(site_index);
    } 


} // namespace Lattice
