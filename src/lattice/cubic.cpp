#include "lattice/cubic.h"


namespace Lattice {

    // high symmetry points in the reciprocal lattice
    const LatticeInt Cubic::GammaPointIndex() const { return this->m_gamma_point_index; }
    const LatticeInt Cubic::XPointIndex() const { return this->m_x_point_index; }
    const LatticeInt Cubic::MPointIndex() const { return this->m_m_point_index; }
    const LatticeInt Cubic::RPointIndex() const { return this->m_r_point_index; }
    const LatticeIntVec& Cubic::DeltaLineIndex() const { return this->m_delta_line_index; }
    const LatticeIntVec& Cubic::ZLineIndex() const { return this->m_z_line_index; }
    const LatticeIntVec& Cubic::SigmaLineIndex() const { return this->m_sigma_line_index; }
    const LatticeIntVec& Cubic::LambdaLineIndex() const { return this->m_lambda_line_index; }
    const LatticeIntVec& Cubic::SLineIndex() const { return this->m_s_line_index; }
    const LatticeIntVec& Cubic::TLineIndex() const { return this->m_t_line_index; }
    
    
    void Cubic::set_lattice_params( const LatticeIntVec& side_length_vec )
    {
        // lattice in three dimension
        assert( (int)side_length_vec.size() == 3 );
        // for cubic lattice, the length of each side should be equal to each other
        assert(    ( side_length_vec[0] == side_length_vec[1] ) 
                && ( side_length_vec[0] == side_length_vec[2] ) );
        assert( side_length_vec[0] >= 2 );

        this->m_space_dim = 3;
        this->m_coordination_number = 6;
        this->m_side_length = side_length_vec[0];
        this->m_space_size = side_length_vec[0] * side_length_vec[1] * side_length_vec[2];
    }


    void Cubic::initial_index2site_table()
    {
        this->m_index2site_table.resize(this->m_space_size, this->m_space_dim);
        for (auto index = 0; index < this->m_space_size; ++index) {
            // map the site index to the site vector (x,y,z)
            this->m_index2site_table(index, 0) = index % this->m_side_length;
            this->m_index2site_table(index, 1) = ( index % (this->m_side_length * this->m_side_length) ) / this->m_side_length;
            this->m_index2site_table(index, 2) = index / (this->m_side_length * this->m_side_length);
        }
    }


    void Cubic::initial_index2momentum_table()
    {
        // the k stars (inequivalent momentum points) in 3d cubic lattice
        // locate in the trirectangular tetrahedron determined by the four vertices 
        // (0,0,0), (pi,0,0), (pi,pi,0) and (pi,pi,pi).
        // note that the point group of 3d cubic lattice is Oh
        const auto half_side_length = std::floor(this->m_side_length/2.0) + 1;
        this->m_num_k_stars = half_side_length * ( half_side_length+1 ) * ( 2*half_side_length+1 ) / 12 
                            + half_side_length * ( half_side_length+1 ) / 4;
        
        // initialize indices of k stars
        this->m_k_stars_index.reserve(this->m_num_k_stars);
        for (auto index = 0; index < this->m_num_k_stars; ++index) {
            this->m_k_stars_index.emplace_back(index);
        }
        
        // initialize index2momentum table
        this->m_index2momentum_table.resize(this->m_num_k_stars, this->m_space_dim);
        int count = 0;
        // index of grids along z axis
        for ( auto k = std::ceil(this->m_side_length/2.0); k <= this->m_side_length; ++k ) {
            // index of grids along x axis
            for ( auto i = k; i <= this->m_side_length; ++i ) {
                // index of grids along y axis
                for ( auto j = k; j <= i; ++j ) {
                    this->m_index2momentum_table.row(count) = Eigen::Vector3d(
                        (double)i/this->m_side_length * 2*M_PI - M_PI, 
                        (double)j/this->m_side_length * 2*M_PI - M_PI,
                        (double)k/this->m_side_length * 2*M_PI - M_PI
                    );
                    count++;
                }
            }
        }
    }


    void Cubic::initial_nearest_neighbour_table()
    {
        // the coordination number for 3d cubic lattice is 6
        // correspondense between the table index and the direction of displacement :
        // 0: (x+1, y, z)    1: (x, y+1, z)    2: (x, y, z+1)
        // 3: (x-1, y, z)    4: (x, y-1, z)    5: (x, y, z-1)
        this->m_nearest_neighbour_table.resize(this->m_space_size, this->m_coordination_number);
        for (int index = 0; index < this->m_space_size; ++index) {
            const auto x = index % this->m_side_length;
            const auto y = ( index % (this->m_side_length * this->m_side_length) ) / this->m_side_length;
            const auto z = index / (this->m_side_length * this->m_side_length);

            this->m_nearest_neighbour_table(index, 0) = ( (x+1)%this->m_side_length ) + this->m_side_length * y 
                                                        + this->m_side_length * this->m_side_length * z;
            this->m_nearest_neighbour_table(index, 3) = ( (x-1+this->m_side_length)%this->m_side_length ) + this->m_side_length * y 
                                                        + this->m_side_length * this->m_side_length * z;
            this->m_nearest_neighbour_table(index, 1) = x + this->m_side_length * ( (y+1)%this->m_side_length ) 
                                                        + this->m_side_length * this->m_side_length * z;
            this->m_nearest_neighbour_table(index, 4) = x + this->m_side_length * ( (y-1+this->m_side_length)%this->m_side_length ) 
                                                        + this->m_side_length * this->m_side_length * z;
            this->m_nearest_neighbour_table(index, 2) = x + this->m_side_length * y
                                                        + this->m_side_length * this->m_side_length * ( (z+1)%this->m_side_length );
            this->m_nearest_neighbour_table(index, 5) = x + this->m_side_length * y
                                                        + this->m_side_length * this->m_side_length * ( (z-1+this->m_side_length)%this->m_side_length );
        }
    }


    void Cubic::initial_displacement_table()
    {
        this->m_displacement_table.resize(this->m_space_size, this->m_space_size);
        for (auto i = 0; i < this->m_space_size; ++i) {
            const auto xi = i % this->m_side_length;
            const auto yi = ( i % (this->m_side_length * this->m_side_length) ) / this->m_side_length;
            const auto zi = i / (this->m_side_length * this->m_side_length);
            
            for (auto j = 0; j < this->m_space_size; ++j) {
                const auto xj = j % this->m_side_length;
                const auto yj = ( j % (this->m_side_length * this->m_side_length) ) / this->m_side_length;
                const auto zj = j / (this->m_side_length * this->m_side_length);

                // displacement pointing from site i to site j
                const auto dx = (xj - xi + this->m_side_length) % this->m_side_length;
                const auto dy = (yj - yi + this->m_side_length) % this->m_side_length;
                const auto dz = (zj - zi + this->m_side_length) % this->m_side_length;
                this->m_displacement_table(i, j) = dx + dy * this->m_side_length + dz * this->m_side_length * this->m_side_length;
            }
        }
    }


    void Cubic::initial_symmetry_points() 
    {
        // high symmetry points of 3d cubic lattice
        // Gamma point:  ( 0,  0,  0)
        // X point:      (pi,  0,  0)
        // M point:      (pi, pi,  0)
        // R point:      (pi, pi, pi)
        this->m_gamma_point_index = 0;
        this->m_m_point_index = ( std::floor(this->m_side_length/2.0) + 1 ) 
                              * ( std::floor(this->m_side_length/2.0) + 2 ) / 2 - 1;
        this->m_x_point_index = this->m_m_point_index - std::floor(this->m_side_length/2.0);
        this->m_r_point_index = this->m_num_k_stars - 1;

        // high symmetry lines of 3d cubic lattice
        // Delta line:   (0,0,0)   ->  (pi,0,0)
        // Z line:       (pi,0,0)  ->  (pi,pi,0)
        // Sigma line:   (0,0,0)   ->  (pi,pi,0)
        // Lambda line:  (0,0,0)   ->  (pi,pi,pi)
        // S line:       (pi,0,0)  ->  (pi,pi,pi)
        // T line:       (pi,pi,0) ->  (pi,pi,pi)
        for ( auto i = 0; i < std::floor(this->m_side_length/2.0)+1; ++i ) {
            this->m_delta_line_index.emplace_back(i*(i+1)/2);
            this->m_z_line_index.emplace_back(this->m_x_point_index+i);
            this->m_sigma_line_index.emplace_back(i*(i+3)/2);

            if ( i == 0 ) { 
                this->m_lambda_line_index.emplace_back(this->m_gamma_point_index);
                this->m_s_line_index.emplace_back(this->m_x_point_index);
                this->m_t_line_index.emplace_back(this->m_m_point_index);
            }
            else {
                this->m_lambda_line_index.emplace_back( 
                        this->m_lambda_line_index.back() 
                        + ( std::floor(this->m_side_length/2.0)+2-i ) 
                        * ( std::floor(this->m_side_length/2.0)+3-i ) / 2 
                    );
                this->m_s_line_index.emplace_back(                         
                        this->m_s_line_index.back() 
                        + ( std::floor(this->m_side_length/2.0)+1-i ) 
                        * ( std::floor(this->m_side_length/2.0)+2-i ) / 2 + 1
                    );
                this->m_t_line_index.emplace_back( 
                        this->m_t_line_index.back() 
                        + ( std::floor(this->m_side_length/2.0)+1-i ) 
                        * ( std::floor(this->m_side_length/2.0)+2-i ) / 2
                    );
            }
        }
    }


    void Cubic::initial_fourier_factor_table()
    {
        // Re( exp(-ikx) ) for lattice site x and momentum k 
        this->m_fourier_factor_table.resize(this->m_space_size, this->m_num_k_stars);
        for (auto x_index = 0; x_index < this->m_space_size; ++x_index) {
            for (auto k_index = 0; k_index < this->m_num_k_stars; ++k_index) {
                // this defines the inner product of a site vector x and a momemtum vector k 
                this->m_fourier_factor_table(x_index, k_index) = cos( 
                        ( - this->m_index2site_table(x_index,0) * this->m_index2momentum_table(k_index,0)
                          - this->m_index2site_table(x_index,1) * this->m_index2momentum_table(k_index,1)
                          - this->m_index2site_table(x_index,2) * this->m_index2momentum_table(k_index,2) )
                );
            }
        }
    }


    void Cubic::initial_hopping_matrix()
    {
        this->m_hopping_matrix.resize(this->m_space_size, this->m_space_size);
        for (auto index = 0; index < this->m_space_size; ++index) {
            // direction 0 for x+1, 1 for y+1 and 2 for z+1
            const int index_xplus1 = this->NearestNeighbour(index, 0);
            const int index_yplus1 = this->NearestNeighbour(index, 1);
            const int index_zplus1 = this->NearestNeighbour(index, 2);

            this->m_hopping_matrix(index, index_xplus1) += 1.0;
            this->m_hopping_matrix(index_xplus1, index) += 1.0;
            this->m_hopping_matrix(index, index_yplus1) += 1.0;
            this->m_hopping_matrix(index_yplus1, index) += 1.0;
            this->m_hopping_matrix(index, index_zplus1) += 1.0;
            this->m_hopping_matrix(index_zplus1, index) += 1.0;
        }
    }


    void Cubic::initial()
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
