#include "lattice/square.h"


namespace Lattice {

    void Square::set_lattice_params(const SideLengthVec& side_length_vec) 
    {   
        // lattice in two dimension
        assert( (int)side_length_vec.size() == 2 );
        // for square lattice, the length of each side should be equal to each other
        assert( side_length_vec[0] == side_length_vec[1] );
        assert( side_length_vec[0] >= 2 );

        this->m_space_dim = 2;
        this->m_coordination_number = 4;
        this->m_side_length = side_length_vec[0];
        this->m_space_size = side_length_vec[0] * side_length_vec[1];
    }


    void Square::initial_index2site_table()
    {
        this->m_index2site_table.resize(this->m_space_size, this->m_space_dim);
        for (int index = 0; index < this->m_space_size; ++index) {
            // map the site index to the site vector (x,y)
            this->m_index2site_table(index, 0) = index % this->m_side_length;
            this->m_index2site_table(index, 1) = index / this->m_side_length;
        }
    }


    void Square::initial_nearest_neighbour_table()
    {   
        // the coordination number for 2d square lattice is 4
        // correspondense between the table index and the direction of displacement :
        // 0: (x+1, y)    1: (x, y+1)
        // 2: (x-1, y)    3: (x, y-1)
        this->m_nearest_neighbour_table.resize(this->m_space_size, this->m_coordination_number);
        for (int index = 0; index < this->m_space_size; ++index) {
            const auto x = index % this->m_side_length;
            const auto y = index / this->m_side_length;

            this->m_nearest_neighbour_table(index, 0) = ((x+1)%this->m_side_length) + this->m_side_length * y;
            this->m_nearest_neighbour_table(index, 2) = ((x-1)%this->m_side_length) + this->m_side_length * y;
            this->m_nearest_neighbour_table(index, 1) = (x%this->m_side_length) + this->m_side_length * ((y+1)%this->m_side_length);
            this->m_nearest_neighbour_table(index, 3) = (x%this->m_side_length) + this->m_side_length * ((y-1)%this->m_side_length);
        }
    }
    

    void Square::initial_hopping_matrix()
    {
        this->m_hopping_matrix.resize(this->m_space_size, this->m_space_size);
        for (auto index = 0; index < this->m_space_size; ++index) {
            // direction 0 for x+1 and 1 for y+1 
            const int index_xplus1 = this->NearestNeighbour(index, 0);
            const int index_yplus1 = this->NearestNeighbour(index, 1);

            this->m_hopping_matrix(index, index_xplus1) += 1.0;
            this->m_hopping_matrix(index_xplus1, index) += 1.0;
            this->m_hopping_matrix(index, index_yplus1) += 1.0;
            this->m_hopping_matrix(index_yplus1, index) += 1.0;
        }
    }


    void Square::initial()
    {
        this->initial_nearest_neighbour_table();
        this->initial_index2site_table();
        this->initial_hopping_matrix();
    }


    // const double Square::product(const std::array<double,2>& vecr, const std::array<double,2>& vecp) {
    //     const auto& [rx, ry] = vecr;
    //     const auto& [px, py] = vecp;
    //     // for square lattice, the two basis are orthogonal
    //     return rx * px + ry * py;
    // }

} // namespace Lattice