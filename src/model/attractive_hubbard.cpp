#include "model/attractive_hubbard.h"
#include "lattice/lattice_base.h"
#include "dqmc_walker.h"
#include "random.h"

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>


namespace Model {
    
    using RealScalar = double;
    using SpaceTimeMat = Eigen::MatrixXd;
    using SpaceSpaceMat = Eigen::MatrixXd;


    const RealScalar AttractiveHubbard::HoppingT() const {
        return this->m_hopping_t;
    }
    
    const RealScalar AttractiveHubbard::ChemicalPotential() const {
        return this->m_chemical_potential;
    }

    const RealScalar AttractiveHubbard::OnSiteU()  const {
        return this->m_onsite_u;
    }


    void AttractiveHubbard::set_model_params (  RealScalar hopping_t, 
                                                RealScalar onsite_u, 
                                                RealScalar chemical_potential  ) 
    {
        assert( hopping_t >= 0.0 );
        assert( onsite_u >= 0.0 );              // abs of the onsite attractive interaction
        this->m_hopping_t = hopping_t;
        this->m_onsite_u = onsite_u;
        this->m_chemical_potential = chemical_potential;
    }


    void AttractiveHubbard::initial_params( const LatticeBase& lattice, const Walker& walker )
    {
        this->m_space_size = lattice.SpaceSize();
        this->m_time_size  = walker.TimeSize();
        const RealScalar time_interval = walker.TimeInterval();

        this->m_alpha = acosh( exp(0.5 * time_interval * this->m_onsite_u) );

        // allocate memory for bosonic fields
        this->m_bosonic_field.resize(this->m_time_size, this->m_space_size);
    }


    void AttractiveHubbard::initial_KV_matrices( const LatticeBase& lattice, const Walker& walker ) 
    {   
        const int space_size = lattice.SpaceSize();
        const RealScalar time_interval = walker.TimeInterval();
        const SpaceSpaceMat chemical_potential_mat = this->m_chemical_potential * SpaceSpaceMat::Identity(space_size,space_size);
        const SpaceSpaceMat Kmat = -this->m_hopping_t * lattice.HoppingMatrix() + chemical_potential_mat;
        
        this->m_expK_mat       = ( -time_interval * Kmat ).exp();
        this->m_inv_expK_mat   = ( +time_interval * Kmat ).exp();
        
        // in general K matrix is symmetrical
        this->m_trans_expK_mat = this->m_expK_mat.transpose();

        // since V is diagonalized in the Hubbard model
        // there is no need to explicitly compute expV
    }


    void AttractiveHubbard::initial( const LatticeBase& lattice, const Walker& walker ) 
    {   
        // initialize model params and allocate memory for bosonic fields
        this->initial_params(lattice, walker);

        // initialize K matrices
        // no need to initialize V matrices because in our model V is diagonalized.
        this->initial_KV_matrices(lattice, walker);
    }
    
    
    void AttractiveHubbard::set_bosonic_fields_to_random() 
    {
        // set configurations of the bosonic fields to random
        const auto time_size = this->m_bosonic_field.rows();
        const auto space_size = this->m_bosonic_field.cols();

        std::bernoulli_distribution bernoulli_dist(0.5);
        for (auto t = 0; t < time_size; ++t) {
            for (auto i = 0; i < space_size; ++i) {
                // for Z2 bosonic field, simply set 1.0 or -1.0
                this->m_bosonic_field(t, i) = bernoulli_dist(Utils::Random::Engine)? +1.0 : -1.0;
            }
        }
    }


    void AttractiveHubbard::update_bosonic_field( TimeIndex time_index, SpaceIndex space_index )
    {
        assert( time_index >= 0 && time_index < this->m_time_size );
        assert( space_index >= 0 && space_index < this->m_space_size );

        // for Z2 bosonic fields, a local update is presented by a local Z2 flip
        this->m_bosonic_field(time_index, space_index) = - this->m_bosonic_field(time_index, space_index);
    }


    const double AttractiveHubbard::get_update_ratio( Walker& walker, TimeIndex time_index, SpaceIndex space_index ) const
    {
        assert( time_index >= 0 && time_index < this->m_time_size );
        assert( space_index >= 0 && space_index < this->m_space_size );

        const Eigen::MatrixXd& green_tt_up = walker.GreenttUp();
        const Eigen::MatrixXd& green_tt_dn = walker.GreenttDn();

        return  exp( 2 * this->m_alpha * this->m_bosonic_field(time_index, space_index) ) 
                * (  1 + ( 1 - green_tt_up(space_index, space_index) ) 
                       * ( exp( -2 * this->m_alpha * this->m_bosonic_field(time_index, space_index) ) - 1 )  )
                * (  1 + ( 1 - green_tt_dn(space_index, space_index) ) 
                       * ( exp( -2 * this->m_alpha * this->m_bosonic_field(time_index, space_index) ) - 1 )  );
    }


    void AttractiveHubbard::update_greens_function( Walker& walker, TimeIndex time_index, SpaceIndex space_index )
    {
        // update the equal-time greens functions 
        // as a consequence of a local Z2 flip of the bosonic fields at (time_index, space_index)
        assert( time_index >= 0 && time_index < this->m_time_size );
        assert( space_index >= 0 && space_index < this->m_space_size );

        Eigen::MatrixXd& green_tt_up = walker.GreenttUp();
        Eigen::MatrixXd& green_tt_dn = walker.GreenttDn();

        // reference:
        //   Quantum Monte Carlo Methods (Algorithms for Lattice Models) Determinant method
        // here we use the sparseness of the matrix \delta
        const double factor_up 
            = ( exp( -2 * this->m_alpha * this->m_bosonic_field(time_index, space_index) ) - 1 )
            / ( 1 + ( 1 - green_tt_up(space_index, space_index) ) 
            * ( exp( -2 * this->m_alpha * this->m_bosonic_field(time_index, space_index) ) - 1 ) );
       
        // for attractive hubbard model, because the spin-up and spin-down parts are coupled 
        // to the bosonic fields in the same way, the model possesses the spin degeneracy, 
        // and in princile it's sufficient to only simulate one specific spin state.
        const double factor_dn = factor_up;
        
        green_tt_up
            -= factor_up * green_tt_up.col(space_index) 
            * ( Eigen::VectorXd::Unit(this->m_space_size, space_index).transpose() - green_tt_up.row(space_index) );
        green_tt_dn
            -= factor_dn * green_tt_dn.col(space_index)
            * ( Eigen::VectorXd::Unit(this->m_space_size, space_index).transpose() - green_tt_dn.row(space_index) );
    }


    void AttractiveHubbard::mult_B_from_left( GreensFunc& green, TimeIndex time_index, Spin spin ) const
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
        this->m_mult_expK_from_left( green );
        for (auto i = 0; i < this->m_space_size; ++i) {
            green.row(i) *= exp( +1.0 * this->m_alpha * this->m_bosonic_field(eff_time_index, i) );
        }
    }


    void AttractiveHubbard::mult_B_from_right( GreensFunc& green, TimeIndex time_index, Spin spin ) const
    {
        // Multiply a dense matrix, specifically a greens function, from the right by B(t)
        //      G  ->  G * B(t) = G * exp( -dt V_sigma(t) ) * exp( -dt K )
        // Matrix G is changed in place.
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        assert( time_index >= 0 && time_index <= this->m_time_size );
        assert( abs(spin) == 1.0 );

        const int eff_time_index = ( time_index == 0 )? this->m_time_size-1 : time_index-1;
        for (auto i = 0; i < this->m_space_size; ++i) {
            green.col(i) *= exp( +1.0 * this->m_alpha * this->m_bosonic_field(eff_time_index, i) );
        }
        this->m_mult_expK_from_right( green );
    }


    void AttractiveHubbard::mult_invB_from_left( GreensFunc& green, TimeIndex time_index, Spin spin ) const
    {
        // Multiply a dense matrix, specifically a greens function, from the left by B(t)^-1
        //      G  ->  B(t)^-1 * G = exp( +dt K ) * exp( +dt V_sigma(t) ) * G
        // Matrix G is changed in place.
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        assert( time_index >= 0 && time_index <= this->m_time_size );
        assert( abs(spin) == 1.0 );

        const int eff_time_index = ( time_index == 0 )? this->m_time_size-1 : time_index-1;
        for (auto i = 0; i < this->m_space_size; ++i) {
            green.row(i) *= exp( -1.0 * this->m_alpha * this->m_bosonic_field(eff_time_index, i) );
        }
        this->m_mult_inv_expK_from_left( green );
    }


    void AttractiveHubbard::mult_invB_from_right( GreensFunc& green, TimeIndex time_index, Spin spin ) const
    {
        // Multiply a dense matrix, specifically a greens function, from the right by B(t)^-1
        //      G  ->  G * B(t)^-1 = G * exp( +dt K ) * exp( +dt V_sigma(t) )
        // Matrix G is changed in place.
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        assert( time_index >= 0 && time_index <= this->m_time_size );
        assert( abs(spin) == 1.0 );

        const int eff_time_index = ( time_index == 0 )? this->m_time_size-1 : time_index-1;
        this->m_mult_inv_expK_from_right( green );
        for (auto i = 0; i < this->m_space_size; ++i) {
            green.col(i) *= exp( -1.0 * this->m_alpha * this->m_bosonic_field(eff_time_index, i) );
        }
    }

    
    void AttractiveHubbard::mult_transB_from_left( GreensFunc& green, TimeIndex time_index, Spin spin ) const
    {
        // Multiply a dense matrix, specifically a greens function, from the left by B(t)^T
        //      G  ->  B(t)^T * G = exp( -dt K )^T * exp( -dt V_sigma(t) ) * G
        // Matrix G is changed in place.
        assert( green.rows() == this->m_space_size && green.cols() == this->m_space_size );
        assert( time_index >= 0 && time_index <= this->m_time_size );
        assert( abs(spin) == 1.0 );

        const int eff_time_index = ( time_index == 0 )? this->m_time_size-1 : time_index-1;
        for (auto i = 0; i < this->m_space_size; ++i) {
            green.row(i) *= exp( +1.0 * this->m_alpha * this->m_bosonic_field(eff_time_index, i) );
        }
        this->m_mult_trans_expK_from_left( green );
    }


} // namespace Model 
