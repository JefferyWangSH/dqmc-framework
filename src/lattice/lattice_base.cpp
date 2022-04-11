#include "lattice/lattice_base.h"

namespace Lattice {

    LatticeBase::LatticeBase(int space_size) {
        this->set_space_size(space_size);
    }

    void LatticeBase::set_space_size(int space_size) { 
        this->m_space_size = space_size; 
        this->m_total_site_num = (int)std::pow(this->m_space_size, this->m_space_dim);
    }
    
    int LatticeBase::SpaceDim() { return this->m_space_dim; }

    int LatticeBase::SpaceSize() { return this->m_space_size; }

    int LatticeBase::TotalSiteNum() { return this->m_total_site_num; }
    
    const Eigen::MatrixXd& LatticeBase::HoppingMatrix() { return this->m_hopping_matrix; }
    
    int LatticeBase::site2index(const std::array<int,2>& site) {
        const auto& [x, y] = site;
        return (x % this->m_space_size) + this->m_space_size * (y % this->m_space_size);
    }

    const std::array<int,2> LatticeBase::index2site(int index) {
        assert( index >= 0 && index < this->m_total_site_num );
        const int x = index % this->m_space_size;
        const int y = index / this->m_space_size;
        return std::array<int,2>{x, y};
    }

    void LatticeBase::initial() {
        this->m_hopping_matrix.resize(this->m_total_site_num, this->m_total_site_num);
        for(int x = 0; x < this->m_space_size; ++x) {
            for(int y = 0; y < this->m_space_size; ++y) {
                const int index = this->site2index({x, y});
                const int index_xplus1 = this->site2index({x+1, y});
                const int index_yplus1 = this->site2index({x, y+1});
                this->m_hopping_matrix(index, index_xplus1) = -1.0;
                this->m_hopping_matrix(index_xplus1, index) = -1.0;
                this->m_hopping_matrix(index, index_yplus1) = -1.0;
                this->m_hopping_matrix(index_yplus1, index) = -1.0;
            }
        }
    }

} // namespace Lattice
