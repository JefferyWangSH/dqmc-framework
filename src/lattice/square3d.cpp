#include "lattice/square3d.h"

namespace Lattice {

    Square3d::Square3d(int space_size) { 
        this->m_space_dim = 3;
        this->m_space_size = space_size; 
        this->m_total_site_num = (int)std::pow(this->m_space_size, this->m_space_dim);  
    }

    void Square3d::initial() {
        this->m_hopping_matrix.resize(this->m_total_site_num, this->m_total_site_num);
        for(int x = 0; x < this->m_space_size; ++x) {
            for(int y = 0; y < this->m_space_size; ++y) {
                for (int z = 0; z < this->m_space_size; ++z) {
                    const int index = this->site2index({x, y, z});
                    const int index_xplus1 = this->site2index({x+1, y, z});
                    const int index_yplus1 = this->site2index({x, y+1, z});
                    const int index_zplus1 = this->site2index({x, y, z+1});
                    this->m_hopping_matrix(index, index_xplus1) = -1.0;
                    this->m_hopping_matrix(index_xplus1, index) = -1.0;
                    this->m_hopping_matrix(index, index_yplus1) = -1.0;
                    this->m_hopping_matrix(index_yplus1, index) = -1.0;
                    this->m_hopping_matrix(index, index_zplus1) = -1.0;
                    this->m_hopping_matrix(index_zplus1, index) = -1.0;
                }
            }
        }
    }

} // namespace Lattice
