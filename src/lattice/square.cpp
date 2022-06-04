#include "lattice/square.h"


namespace Lattice {

    // high symmetry points in the reciprocal lattice
    const LatticeInt Square::GammaPointIndex() const { return this->m_gamma_point_index; }
    const LatticeInt Square::XPointIndex() const { return this->m_x_point_index; }
    const LatticeInt Square::MPointIndex() const { return this->m_m_point_index; }
    const LatticeIntVec& Square::DeltaLineIndex() const { return this->m_delta_line_index; }
    const LatticeIntVec& Square::ZLineIndex() const { return this->m_z_line_index; }
    const LatticeIntVec& Square::SigmaLineIndex() const { return this->m_sigma_line_index; }
    const LatticeIntVec& Square::Gamma2X2M2GammaLoopIndex() const { return this->m_gamma2x2m2gamma_loop_index; }


    void Square::set_lattice_params(const LatticeIntVec& side_length_vec) 
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
        for (auto index = 0; index < this->m_space_size; ++index) {
            // map the site index to the site vector (x,y)
            this->m_index2site_table(index, 0) = index % this->m_side_length;
            this->m_index2site_table(index, 1) = index / this->m_side_length;
        }
    }

    
    void Square::initial_index2momentum_table()
    {
        // k stars (inequivalent momentum points) in 2d square lattice
        // locate in the zone surrounded by loop (0,0) -> (pi,0) -> (pi,pi) -> (0,0).
        // note that the point group of 2d sqaure lattice is C4v
        this->m_num_k_stars = (std::floor(this->m_side_length/2.0)+1)*(std::floor(this->m_side_length/2.0)+2)/2;
        
        // initialize indices of k stars
        this->m_k_stars_index.reserve(this->m_num_k_stars);
        for (auto index = 0; index < this->m_num_k_stars; ++index) {
            this->m_k_stars_index.emplace_back(index);
        }
        
        // initialize index2momentum table
        this->m_index2momentum_table.resize(this->m_num_k_stars, this->m_space_dim);
        int count = 0; 
        for (auto i = std::ceil(this->m_side_length/2.0); i <= this->m_side_length; ++i) {
            for (auto j = std::ceil(this->m_side_length/2.0); j <= i; ++j) {
                this->m_index2momentum_table.row(count) = Eigen::Vector2d(  
                    (double)i/this->m_side_length * 2*M_PI - M_PI, (double)j/this->m_side_length * 2*M_PI - M_PI
                );
                count++;
            }
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
            this->m_nearest_neighbour_table(index, 1) = x + this->m_side_length * ((y+1)%this->m_side_length);
            this->m_nearest_neighbour_table(index, 3) = x + this->m_side_length * ((y-1)%this->m_side_length);
        }
    }


    void Square::initial_displacement_table()
    {
        this->m_displacement_table.resize(this->m_space_size, this->m_space_size);
        for (auto i = 0; i < this->m_space_size; ++i) {
            const auto xi = i % this->m_side_length;
            const auto yi = i / this->m_side_length;
            
            for (auto j = 0; j < this->m_space_size; ++j) {
                const auto xj = j % this->m_side_length;
                const auto yj = j / this->m_side_length;

                // displacement pointing from site i to site j
                const auto dx = (xj - xi + this->m_side_length) % this->m_side_length;
                const auto dy = (yj - yi + this->m_side_length) % this->m_side_length;
                this->m_displacement_table(i, j) = dx + dy * this->m_side_length;
            }
        }
    }
            

    void Square::initial_symmetry_points() 
    {
        // high symmetry points of 2d square lattice
        // Gamma point:  (0,  0)
        // X point:      (pi, 0)
        // M point:      (pi, pi)
        this->m_gamma_point_index = 0;
        this->m_x_point_index     = this->m_num_k_stars - std::floor(this->m_side_length/2.0) - 1;
        this->m_m_point_index     = this->m_num_k_stars - 1; 

        // high symmetry lines of 2d square lattice
        // Delta line:   (0,0)  ->  (pi,0)
        // Z line:       (pi,0) ->  (pi,pi)
        // Sigma line:   (0,0)  ->  (pi,pi)
        this->m_delta_line_index.reserve(std::floor(this->m_side_length/2.0)+1);
        this->m_z_line_index.reserve(std::floor(this->m_side_length/2.0)+1);
        this->m_sigma_line_index.reserve(std::floor(this->m_side_length/2.0)+1);
        for (auto i = 0; i < std::floor(this->m_side_length/2.0)+1; ++i) {
            this->m_delta_line_index.emplace_back(i*(i+1)/2);
            this->m_z_line_index.emplace_back(this->m_x_point_index+i);
            this->m_sigma_line_index.emplace_back(i*(i+3)/2);
        }

        // loop: (0,0) -> (pi,0) -> (pi,pi) -> (0,0)
        this->m_gamma2x2m2gamma_loop_index.reserve(3*(this->m_side_length-std::ceil(this->m_side_length/2.0)));
        for (auto i = 0; i < std::floor(this->m_side_length/2.0); ++i) {
            // along (0,0) -> (pi,0) direation
            this->m_gamma2x2m2gamma_loop_index.emplace_back(i*(i+1)/2);
        }
        for (auto i = 0; i < std::floor(this->m_side_length/2.0); ++i) {
            // along (pi,0) -> (pi,pi) direction
            this->m_gamma2x2m2gamma_loop_index.emplace_back(this->m_x_point_index+i);
        }
        for (auto i = std::floor(this->m_side_length/2.0); i >= 1; --i) {
            // along (pi,pi) -> (0,0) direction
            this->m_gamma2x2m2gamma_loop_index.emplace_back(i*(i+3)/2);
        }
    }


    void Square::initial_fourier_factor_table()
    {
        // Re( exp(-ikx) ) for lattice site x and momentum k 
        this->m_fourier_factor_table.resize(this->m_space_size, this->m_num_k_stars);
        for (auto x_index = 0; x_index < this->m_space_size; ++x_index) {
            for (auto k_index = 0; k_index < this->m_num_k_stars; ++k_index) {
                // this defines the inner product of a site vector x and a momemtum vector k 
                this->m_fourier_factor_table(x_index, k_index) = cos( 
                        ( - this->m_index2site_table(x_index,0) * this->m_index2momentum_table(k_index,0)
                          - this->m_index2site_table(x_index,1) * this->m_index2momentum_table(k_index,1) )
                );
            }
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
        // avoid multiple initialization
        if ( !this->m_initial_status ) {
            this->initial_index2site_table();
            this->initial_index2momentum_table();

            this->initial_nearest_neighbour_table();
            this->initial_displacement_table();
            this->initial_symmetry_points();
            this->initial_fourier_factor_table();

            this->initial_hopping_matrix();  

            this->m_initial_status = true;
        }
    }


} // namespace Lattice