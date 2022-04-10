#include "lattice/square2d.h"

namespace Lattice {

    Square2d::Square2d(int space_size) { 
        this->m_space_dim = 2;
        this->m_space_size = space_size;
        this->m_total_site_num = (int)std::pow(this->m_space_size, this->m_space_dim);  
    }

    void Square2d::initial() {
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
