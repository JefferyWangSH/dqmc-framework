#include "model/repulsive_hubbard_cbsquare2d.h"
#include "lattice/square2d.h"
#include "dqmc_walker.h"


namespace Model {

    void RepulsiveHubbardCbSquare2d::initial_checkerboard(const LatticeBase& lattice, const Walker& walker) 
    {   
        // ensure the lattice is 2d square type
        assert( dynamic_cast<const Lattice::Square2d*>(&lattice) != nullptr );

        const int side_length = lattice.SpaceSize();
        const int space_size = lattice.TotalSiteNum();
        const RealScalar time_interval = walker.TimeInterval();

        this->m_checkerboard.set_params(side_length, space_size, time_interval, this->m_hopping_t, this->m_chemical_potential);
        this->m_checkerboard.initial();
    }


    void RepulsiveHubbardCbSquare2d::initial(const LatticeBase& lattice, const Walker& walker) 
    {   
        // initialize parmas of base class
        RepulsiveHubbard::initial_params(lattice, walker);

        // initialize checkerboard class
        this->initial_checkerboard(lattice, walker);
    }


    void RepulsiveHubbardCbSquare2d::mult_B_from_left( GreensFunc& green, TimeIndex time_index, Spin spin ) const
    {
        // Multiply a dense matrix, specifically a greens function, from the left by B(t)
        //      G  ->  B(t) * G = exp( -dt V_sigma(t) ) * exp( -dt K ) * G
        // Matrix G is changed in place.
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        assert( time_index >= 0 && time_index <= this->m_time_size );
        // 1.0 for spin up and -1.0 for spin down
        assert( abs(spin) == 1.0 );

        // due to the periodical boundary condition (PBC)
        // the time slice labeled by 0 actually corresponds to slice tau = beta
        const int eff_time_index = ( time_index == 0 )? this->m_time_size-1 : time_index-1;        
        this->m_checkerboard.mult_expK_from_left(green);
        for (auto i = 0; i < this->m_space_size; ++i) {
            green.row(i) *= exp( +spin * this->m_alpha * this->m_bosonic_field(eff_time_index, i) );
        }
    }


    void RepulsiveHubbardCbSquare2d::mult_B_from_right( GreensFunc& green, TimeIndex time_index, Spin spin ) const
    {
        // Multiply a dense matrix, specifically a greens function, from the right by B(t)
        //      G  ->  G * B(t) = G * exp( -dt V_sigma(t) ) * exp( -dt K )
        // Matrix G is changed in place.
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        assert( time_index >= 0 && time_index <= this->m_time_size );
        assert( abs(spin) == 1.0 );

        const int eff_time_index = ( time_index == 0 )? this->m_time_size-1 : time_index-1;
        for (auto i = 0; i < this->m_space_size; ++i) {
            green.col(i) *= exp( +spin * this->m_alpha * this->m_bosonic_field(eff_time_index, i) );
        }
        this->m_checkerboard.mult_expK_from_right(green);
    }


    void RepulsiveHubbardCbSquare2d::mult_invB_from_left( GreensFunc& green, TimeIndex time_index, Spin spin ) const
    {
        // Multiply a dense matrix, specifically a greens function, from the left by B(t)^-1
        //      G  ->  B(t)^-1 * G = exp( +dt K ) * exp( +dt V_sigma(t) ) * G
        // Matrix G is changed in place.
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        assert( time_index >= 0 && time_index <= this->m_time_size );
        assert( abs(spin) == 1.0 );

        const int eff_time_index = ( time_index == 0 )? this->m_time_size-1 : time_index-1;
        for (auto i = 0; i < this->m_space_size; ++i) {
            green.row(i) *= exp( -spin * this->m_alpha * this->m_bosonic_field(eff_time_index, i) );
        }
        this->m_checkerboard.mult_inv_expK_from_left(green);
    }


    void RepulsiveHubbardCbSquare2d::mult_invB_from_right( GreensFunc& green, TimeIndex time_index, Spin spin ) const
    {
        // Multiply a dense matrix, specifically a greens function, from the right by B(t)^-1
        //      G  ->  G * B(t)^-1 = G * exp( +dt K ) * exp( +dt V_sigma(t) )
        // Matrix G is changed in place.
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        assert( time_index >= 0 && time_index <= this->m_time_size );
        assert( abs(spin) == 1.0 );

        const int eff_time_index = ( time_index == 0 )? this->m_time_size-1 : time_index-1;
        this->m_checkerboard.mult_inv_expK_from_right(green);
        for (auto i = 0; i < this->m_space_size; ++i) {
            green.col(i) *= exp( -spin * this->m_alpha * this->m_bosonic_field(eff_time_index, i) );
        }
    }

    
    void RepulsiveHubbardCbSquare2d::mult_transB_from_left( GreensFunc& green, TimeIndex time_index, Spin spin ) const
    {
        // Multiply a dense matrix, specifically a greens function, from the left by B(t)^T
        //      G  ->  B(t)^T * G = exp( -dt K )^T * exp( -dt V_sigma(t) ) * G
        // Matrix G is changed in place.
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        assert( time_index >= 0 && time_index <= this->m_time_size );
        assert( abs(spin) == 1.0 );

        const int eff_time_index = ( time_index == 0 )? this->m_time_size-1 : time_index-1;
        for (auto i = 0; i < this->m_space_size; ++i) {
            green.row(i) *= exp( +spin * this->m_alpha * this->m_bosonic_field(eff_time_index, i) );
        }
        this->m_checkerboard.mult_trans_expK_from_left(green);
    }

} // namespace Model