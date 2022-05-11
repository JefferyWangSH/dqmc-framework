#ifndef LATTICE_BASE_H
#define LATTICE_BASE_H
#pragma once

/**
  *  This header file defines the pure virtual abstract class Lattice::LatticeBase 
  *  for describing the space-discreted lattice where the quantum systems live.
  */


#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


namespace Lattice {

    using LatticeBool = bool;
    using LatticeInt = int;
    using LatticeDouble = double;
    using LatticeIntVec = std::vector<int>;
    using MatrixDouble = Eigen::MatrixXd;
    using VectorDouble = Eigen::VectorXd;
    using MatrixInt = Eigen::MatrixXi;
    using VectorInt = Eigen::VectorXi;


    // -------------------------- Pure virtual base class Lattice::LatticeBase ----------------------------
    class LatticeBase {
        protected:

            LatticeBool m_initial_status{false};        // status of initialization

            LatticeInt  m_space_dim{};                  // dimension of the space
            LatticeInt  m_side_length{};                // side length of the lattice
            LatticeInt  m_space_size{};                 // total number of lattice sites
            LatticeInt  m_coordination_number{};        // coordination number of the lattice
            LatticeInt  m_num_k_stars{};                // number of k stars ( inequivalent momentum points )

            // hopping matrix, depending only on the topology of lattice
            // hopping constants are normalized to 1.0 .
            MatrixDouble m_hopping_matrix{};

            // Matrix structure for storing the nearest neighbours of each lattice site
            // the matrix shape is SpaceSize * Coordination number
            MatrixInt m_nearest_neighbour_table{};
            // todo: next nearest neighbours
            // MatrixInt m_next_nearest_neighbour_table{};
            
            // Matrix structure for storing the map from site index to the site vector
            // with the shape of SpaceSize * SpaceDim
            MatrixInt m_index2site_table{};
            
            // table of the displacement between any two sites i and j, pointing from i to j
            // the displacement is represented by a site index due to the periodic boundary condition,
            // and the shape of the table is SpaceSize * SpaceSize.
            MatrixInt m_displacement_table{};
            
            // the map from momentum index to the lattice momentum in the reciprocal lattice
            // the number of rows should be equal to the number of inequivalent momentum points (k stars),
            // and the columns is the space dimension.
            MatrixDouble m_index2momentum_table{};

            // table of fourier transformation factor
            // e.g. Re( exp(ikx) ) for lattice site x and momentum k 
            MatrixDouble m_fourier_factor_table{};

            // some higher symmetric points in the reciprocal lattice
            // todo: move these to the derived class of specialized lattice
            LatticeInt m_gamma_point_index{};
            LatticeInt m_m_point_index{};
            LatticeInt m_x_point_index{};
            LatticeIntVec m_delta_line_index{};
            LatticeIntVec m_z_line_index{};
            LatticeIntVec m_sigma_line_index{};
            LatticeIntVec m_gamma2x2m2gamma_loop_index{};   // indices of the defined loop
            LatticeIntVec m_k_stars_index{};                // all inequivalent momentum points


        public:
            LatticeBase() = default;

            // ---------------------------- Set up lattice parameters ------------------------------
            // read lattice params from a vector of side lengths
            virtual void set_lattice_params( const LatticeIntVec& side_length_vec ) = 0;

            
            // --------------------------------- Interfaces ----------------------------------------
            
            const LatticeBool InitialStatus()      const ;
            const LatticeInt  SpaceDim()           const ;
            const LatticeInt  SpaceSize()          const ;
            const LatticeInt  SideLength()         const ;
            const LatticeInt  CoordinationNumber() const ;
            const LatticeInt  kStarsNum()          const ;

            // some symmetric points
            const LatticeInt GammaPointIndex()     const ;
            const LatticeInt MPointIndex()         const ;
            const LatticeInt XPointIndex()         const ;
            const LatticeIntVec& DeltaLineIndex()  const ;
            const LatticeIntVec& ZLineIndex()      const ;
            const LatticeIntVec& SigmaLineIndex()  const ;
            const LatticeIntVec& kStarsIndex()     const ;
            const LatticeIntVec& Gamma2X2M2GammaLoopIndex() const ;

            const MatrixDouble& HoppingMatrix() const ;
            const MatrixDouble& FourierFactor() const ;
            const LatticeInt    NearestNeighbour ( const LatticeInt site_index, const LatticeInt direction ) const ;
            const LatticeInt    Displacement     ( const LatticeInt site1_index, const LatticeInt site2_index ) const ;
            const LatticeDouble FourierFactor    ( const LatticeInt site_index, const LatticeInt momentum_index ) const ;

            const VectorInt     Index2Site( const LatticeInt site_index ) const ;
            const LatticeInt    Index2Site( const LatticeInt site_index, const LatticeInt axis ) const ;
            const VectorDouble  Index2Momentum( const LatticeInt momentum_index ) const ;
            const LatticeDouble Index2Momentum( const LatticeInt momentum_index, const LatticeInt axis ) const ;


            // -------------------------------- Initializations ------------------------------------
            
            virtual void initial() = 0;
            virtual void initial_hopping_matrix()          = 0;
            virtual void initial_index2site_table()        = 0;
            virtual void initial_nearest_neighbour_table() = 0;
            virtual void initial_displacement_table()      = 0;
            virtual void initial_index2momentum_table()    = 0;
            virtual void initial_symmetric_points()        = 0;
            virtual void initial_fourier_factor_table()    = 0;

    };

} // namespace Lattice

#endif // LATTICE_BASE_H