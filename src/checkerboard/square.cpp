#include "checkerboard/square.h"
#include "lattice/square.h"
#include "model/model_base.h"
#include "dqmc_walker.h"

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <unsupported/Eigen/MatrixFunctions>


namespace CheckerBoard {

    void Square::set_checkerboard_params( const LatticeBase& lattice, 
                                          const ModelBase& model, 
                                          const DqmcWalker& walker ) 
    {   
        // make sure that the lattice class is of type Lattice::Square
        assert( dynamic_cast<const Lattice::Square*>(&lattice) != nullptr );
        assert( lattice.SpaceSize() >= 2 );
        this->m_side_length = lattice.SideLength();
        this->m_space_size = lattice.SpaceSize();
        this->m_time_interval = walker.TimeInterval();
        this->m_hopping_t = model.HoppingT();
        this->m_chemical_potential = model.ChemicalPotential();
    }


    void Square::initial() 
    {
        // construct the exponent of hopping matrix in a reduced 4*4 Hilbert-space
        Eigen::Matrix4d reduced_hopping_mat;
        double ch, sh;

        ch = cosh(this->m_time_interval * this->m_hopping_t);
        sh = sinh(this->m_time_interval * this->m_hopping_t);
        reduced_hopping_mat <<  ch * ch, ch * sh, ch * sh, sh * sh,
                                ch * sh, ch * ch, sh * sh, ch * sh,
                                ch * sh, sh * sh, ch * ch, ch * sh,
                                sh * sh, ch * sh, ch * sh, ch * ch;
        // factor 0.5 comes from the double counting of sites
        this->m_expK_plaquette = exp( 0.5*this->m_time_interval*this->m_chemical_potential ) * reduced_hopping_mat;

        ch = cosh(-this->m_time_interval * this->m_hopping_t);
        sh = sinh(-this->m_time_interval * this->m_hopping_t);
        reduced_hopping_mat <<  ch * ch, ch * sh, ch * sh, sh * sh,
                                ch * sh, ch * ch, sh * sh, ch * sh,
                                ch * sh, sh * sh, ch * ch, ch * sh,
                                sh * sh, ch * sh, ch * sh, ch * ch;
        this->m_inv_expK_plaquette = exp( -0.5*this->m_time_interval*this->m_chemical_potential ) * reduced_hopping_mat;
    }


    // Ref: Max H. Gerlach, 2017
    // (x,y) counts the upper left conner of a plaquette.
    // plaquette: labeled by site i (x, y), the upper-left conner
    //   i ---- j
    //   |      |
    //   |      |
    //   k ---- l
    // and the corresponding effective hopping matrix (4*4) reads
    //   0.0, 1.0, 1.0, 0.0,
    //   1.0, 0.0, 0.0, 1.0,
    //   1.0, 0.0, 0.0, 1.0,
    //   0.0, 1.0, 1.0, 0.0.

    void Square::mult_expK_plaquette_from_left( Matrix &matrix, const Site& site ) const
    {
        assert( matrix.rows() == this->m_space_size && matrix.cols() == this->m_space_size );
        assert( site.size() == 2 );
        const int x = site[0];
        const int y = site[1];
        const int index_xy          = ( x%this->m_side_length ) + this->m_side_length * ( y%this->m_side_length );
        const int index_xy_right    = ( (x+1)%this->m_side_length ) + this->m_side_length * ( y%this->m_side_length );
        const int index_xy_down     = ( x%this->m_side_length ) + this->m_side_length * ( (y+1)%this->m_side_length );
        const int index_xy_diagonal = ( (x+1)%this->m_side_length ) + this->m_side_length * ( (y+1)%this->m_side_length );
        
        const std::array<int,4> indexes = { index_xy, index_xy_right, index_xy_down, index_xy_diagonal };
        matrix(indexes, Eigen::all) = this->m_expK_plaquette * matrix(indexes, Eigen::all);
    }


    void Square::mult_inv_expK_plaquette_from_left( Matrix &matrix, const Site& site ) const
    {
        assert( matrix.rows() == this->m_space_size && matrix.cols() == this->m_space_size );
        assert( site.size() == 2 );
        const int x = site[0];
        const int y = site[1];
        const int index_xy          = ( x%this->m_side_length ) + this->m_side_length * ( y%this->m_side_length );
        const int index_xy_right    = ( (x+1)%this->m_side_length ) + this->m_side_length * ( y%this->m_side_length );
        const int index_xy_down     = ( x%this->m_side_length ) + this->m_side_length * ( (y+1)%this->m_side_length );
        const int index_xy_diagonal = ( (x+1)%this->m_side_length ) + this->m_side_length * ( (y+1)%this->m_side_length );
        
        const std::array<int,4> indexes = { index_xy, index_xy_right, index_xy_down, index_xy_diagonal };
        matrix(indexes, Eigen::all) = this->m_inv_expK_plaquette * matrix(indexes, Eigen::all);
    }


    void Square::mult_expK_plaquette_from_right( Matrix &matrix, const Site& site ) const
    {
        assert( matrix.rows() == this->m_space_size && matrix.cols() == this->m_space_size );
        assert( site.size() == 2 );
        const int x = site[0];
        const int y = site[1];
        const int index_xy          = ( x%this->m_side_length ) + this->m_side_length * ( y%this->m_side_length );
        const int index_xy_right    = ( (x+1)%this->m_side_length ) + this->m_side_length * ( y%this->m_side_length );
        const int index_xy_down     = ( x%this->m_side_length ) + this->m_side_length * ( (y+1)%this->m_side_length );
        const int index_xy_diagonal = ( (x+1)%this->m_side_length ) + this->m_side_length * ( (y+1)%this->m_side_length );
        
        const std::array<int,4> indexes = { index_xy, index_xy_right, index_xy_down, index_xy_diagonal };
        matrix(Eigen::all, indexes) = matrix(Eigen::all, indexes) * this->m_expK_plaquette;
    }

    
    void Square::mult_inv_expK_plaquette_from_right( Matrix &matrix, const Site& site ) const
    {
        assert( matrix.rows() == this->m_space_size && matrix.cols() == this->m_space_size );
        assert( site.size() == 2 );
        const int x = site[0];
        const int y = site[1];
        const int index_xy          = ( x%this->m_side_length ) + this->m_side_length * ( y%this->m_side_length );
        const int index_xy_right    = ( (x+1)%this->m_side_length ) + this->m_side_length * ( y%this->m_side_length );
        const int index_xy_down     = ( x%this->m_side_length ) + this->m_side_length * ( (y+1)%this->m_side_length );
        const int index_xy_diagonal = ( (x+1)%this->m_side_length ) + this->m_side_length * ( (y+1)%this->m_side_length );
        
        const std::array<int,4> indexes = { index_xy, index_xy_right, index_xy_down, index_xy_diagonal };
        matrix(Eigen::all, indexes) = matrix(Eigen::all, indexes) * this->m_inv_expK_plaquette;
    }


    void Square::mult_expK_from_left( Matrix &matrix ) const
    {
        // checkerboard breakups are only supported for lattices with even side length
        assert( this->m_side_length % 2 == 0 );

        // sublattice B
        for (auto x = 1; x < this->m_side_length; x+=2) {
            for (auto y = 1; y < this->m_side_length; y+=2) {
                this->mult_expK_plaquette_from_left( matrix, {x,y} );
            }
        }
        // sublattice A
        for (int x = 0; x < this->m_side_length; x+=2) {
            for (int y = 0; y < this->m_side_length; y+=2) {
                this->mult_expK_plaquette_from_left( matrix, {x,y} );
            }
        }
    }

    
    void Square::mult_inv_expK_from_left( Matrix &matrix ) const
    {
        // checkerboard breakups are only supported for lattices with even side length
        assert( this->m_side_length % 2 == 0 );

        // sublattice A
        for (int x = 0; x < this->m_side_length; x+=2) {
            for (int y = 0; y < this->m_side_length; y+=2) {
                this->mult_inv_expK_plaquette_from_left( matrix, {x,y} );
            }
        }
        // sublattice B
        for (auto x = 1; x < this->m_side_length; x+=2) {
            for (auto y = 1; y < this->m_side_length; y+=2) {
                this->mult_inv_expK_plaquette_from_left( matrix, {x,y} );
            }
        }
    }

    
    void Square::mult_expK_from_right( Matrix &matrix ) const
    {
        // checkerboard breakups are only supported for lattices with even side length
        assert( this->m_side_length % 2 == 0 );

        // sublattice A
        for (int x = 0; x < this->m_side_length; x+=2) {
            for (int y = 0; y < this->m_side_length; y+=2) {
                this->mult_expK_plaquette_from_right( matrix, {x,y} );
            }
        }
        // sublattice B
        for (auto x = 1; x < this->m_side_length; x+=2) {
            for (auto y = 1; y < this->m_side_length; y+=2) {
                this->mult_expK_plaquette_from_right( matrix, {x,y} );
            }
        }
    }

    
    void Square::mult_inv_expK_from_right( Matrix &matrix ) const
    {
        // checkerboard breakups are only supported for lattices with even side length
        assert( this->m_side_length % 2 == 0 );

        // sublattice B
        for (auto x = 1; x < this->m_side_length; x+=2) {
            for (auto y = 1; y < this->m_side_length; y+=2) {
                this->mult_inv_expK_plaquette_from_right( matrix, {x,y} );
            }
        }
        // sublattice A
        for (int x = 0; x < this->m_side_length; x+=2) {
            for (int y = 0; y < this->m_side_length; y+=2) {
                this->mult_inv_expK_plaquette_from_right( matrix, {x,y} );
            }
        }
    }

    void Square::mult_trans_expK_from_left( Matrix &matrix ) const
    {   
        // checkerboard breakups are only supported for lattices with even side length
        assert( this->m_side_length % 2 == 0 );

        // sublattice A
        for (int x = 0; x < this->m_side_length; x+=2) {
            for (int y = 0; y < this->m_side_length; y+=2) {
                this->mult_expK_plaquette_from_left( matrix, {x,y} );
            }
        }
        // sublattice B
        for (auto x = 1; x < this->m_side_length; x+=2) {
            for (auto y = 1; y < this->m_side_length; y+=2) {
                this->mult_expK_plaquette_from_left( matrix, {x,y} );
            }
        }
    }


} // namespace CheckerBoard
