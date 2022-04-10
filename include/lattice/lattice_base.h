#ifndef LATTICE_BASE_H
#define LATTICE_BASE_H
#pragma once

/**
  *  This header file defines the base class for the space-discreted lattice
  */


#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

namespace Lattice {

    class LatticeBase {
        protected:
            int m_space_dim{};
            int m_space_size{};
            int m_total_site_num{};

            // hopping matrix, depending only on the topology of lattice
            Eigen::MatrixXd m_hopping_matrix{};
            
        public:
            LatticeBase() = default;

            void set_space_size(int space_size);

            int site2index(const std::vector<int>& site);
            const std::vector<int> index2site(int index);

            int SpaceDim();
            int SpaceSize();
            int TotalSiteNum();
            const Eigen::MatrixXd& HoppingMatrix();

            virtual void initial() = 0;
    };

} // namespace Lattice

#endif // LATTICE_BASE_H