#include "lattice/lattice_base.h"

namespace Lattice {

    void LatticeBase::set_space_size(int space_size) { 
        this->m_space_size = space_size; 
        this->m_total_site_num = (int)std::pow(this->m_space_size, this->m_space_dim);
    }

    int LatticeBase::site2index(const std::vector<int>& site) {
        assert( this->m_space_dim == (int)site.size() );
        int index = 0;
        for (int dim = 0; dim < this->m_space_dim; ++dim) {
            index += std::pow(this->m_space_size, dim) * (site[dim]%this->m_space_size);
        }
        return index;
    }

    const std::vector<int> LatticeBase::index2site(int index) {
        assert( index >= 0 && index < this->m_total_site_num );
        std::vector<int> site;
        site.reserve(this->m_space_dim);
        for (int dim = 0; dim < this->m_space_dim; ++dim) {
            site.emplace_back( (index % (int)std::pow(this->m_space_size,dim+1)) / (int)pow(this->m_space_size,dim) );
        }
        return site;
    }

    int LatticeBase::SpaceDim() { return this->m_space_dim; }

    int LatticeBase::SpaceSize() { return this->m_space_size; }

    int LatticeBase::TotalSiteNum() { return this->m_total_site_num; }
    
    const Eigen::MatrixXd& LatticeBase::HoppingMatrix() { return this->m_hopping_matrix; }

} // namespace Lattice
