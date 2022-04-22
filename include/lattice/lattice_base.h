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

    using LatticeInt = int;
    using MatrixDouble = Eigen::MatrixXd;
    using MatrixInt = Eigen::MatrixXi;
    using VectorInt = Eigen::VectorXi;
    using SideLengthVec = std::vector<int>;


    // -------------------------- Pure virtual base class Lattice::LatticeBase ----------------------------
    class LatticeBase {
        protected:
            LatticeInt m_space_dim{};          // dimension of the space
            LatticeInt m_side_length{};        // side length of the lattice
            LatticeInt m_space_size{};         // total number of lattice sites

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
            
            // // todo:
            // // momentum in the reciprocal lattice
            // MatrixInt m_index2momentum_table{};


        public:
            LatticeBase() = default;

            // ---------------------------- Set up lattice parameters ------------------------------
            // read lattice params from a vector of side lengths
            virtual void set_lattice_params(const SideLengthVec& side_length_vec) = 0;

            
            // --------------------------------- Interfaces ----------------------------------------
            
            const LatticeInt SpaceDim()   const ;
            const LatticeInt SpaceSize()  const ;
            const LatticeInt SideLength() const ;

            const MatrixDouble& HoppingMatrix() const ;
            const LatticeInt NearestNeighbour(const LatticeInt site_index, const LatticeInt direction) const ;
            const VectorInt NearestNeighbour(const LatticeInt site_index) const ;
            const VectorInt index2site(const LatticeInt site_index) const ;


            // -------------------------------- Initializations ------------------------------------
            
            virtual void initial() = 0;
            virtual void initial_hopping_matrix()          = 0;
            virtual void initial_index2site_table()        = 0;
            virtual void initial_nearest_neighbour_table() = 0;
            
            // todo
            // virtual void initial_lattice_momentum_table()  = 0;


            // -------------------------- Definition of vector products ----------------------------
            // // vector product of space vector r and moemntum p, depending on the geometry of lattice
            // virtual const double product(const std::array<double,2>& vecr, const std::array<double,2>& vecp) = 0;

    };

} // namespace Lattice

#endif // LATTICE_BASE_H