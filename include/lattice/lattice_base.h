#ifndef LATTICE_BASE_H
#define LATTICE_BASE_H
#pragma once

/**
  *  This header file defines the abstract base class Lattice::LatticeBase 
  *  for describing the space-discreted lattice where the quantum systems live.
  */


#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

namespace Lattice {

        class LatticeBase {
        protected:
            int m_space_dim{2};
            int m_space_size{};
            int m_total_site_num{};

            // hopping matrix, depending only on the topology of lattice
            Eigen::MatrixXd m_hopping_matrix{};
            
        public:
            LatticeBase() = default;
            LatticeBase(int space_size);

            void set_space_size(int space_size);

            int SpaceDim();
            int SpaceSize();
            int TotalSiteNum();
            const Eigen::MatrixXd& HoppingMatrix();

            // todo: maybe replace vector with x and y
            int site2index(const std::array<int,2>& site);
            const std::array<int,2> index2site(int index);

            // initialize class, especially generating hopping matrix
            void initial();

            // vector product of space vector r and moemntum p, depending on the geometry of lattice
            virtual double product(const std::array<double,2>& vecr, const std::array<double,2>& vecp) = 0;

            // // todo: add two space sites according to the topology of lattice, return index of site
            // int add(int index1, int index2);
            
            // // todo: move one step along direction of the 1st space basis, return index of site
            // int move1step(int index);
    };

} // namespace Lattice

#endif // LATTICE_BASE_H