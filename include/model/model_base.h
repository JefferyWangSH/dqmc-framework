#ifndef MODEL_BASE_H
#define MODEL_BASE_H
#pragma once


/**
  *  This header file defines the abstract base class Model::ModelBase, which is pure virtual,
  *  for describing different kinds of quantum models.
  */  

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


// forward declaration
namespace QuantumMonteCarlo {
    class DqmcWalker;
}
namespace Lattice {
    class LatticeBase;
}

namespace Model {

    // -------------------------- Abstract base class Model::ModelBase -----------------------------
    class ModelBase {
        protected:
            using SpaceIndex = int;
            using TimeIndex = int;
            using Spin = int;
            using Walker = QuantumMonteCarlo::DqmcWalker;
            using Lattice = Lattice::LatticeBase;
            using GreensFunc = Eigen::MatrixXd;

        public:
            ModelBase() = default;
        
        protected:
        
            virtual void initial(const Lattice& lattice) = 0;

            virtual double get_update_radio(const Walker& walker, SpaceIndex space_index, TimeIndex time_index) = 0;

            virtual void local_update(const Walker& walker, SpaceIndex space_index, TimeIndex time_index) = 0;

            virtual void mult_B_from_left(GreensFunc& green, const Walker& walker, TimeIndex time_index, Spin spin) = 0;
            virtual void mult_B_from_right(GreensFunc& green, const Walker& walker, TimeIndex time_index, Spin spin) = 0;
            virtual void mult_invB_from_left(GreensFunc& green, const Walker& walker, TimeIndex time_index, Spin spin) = 0;
            virtual void mult_invB_from_right(GreensFunc& green, const Walker& walker, TimeIndex time_index, Spin spin) = 0;
            virtual void mult_transB_from_left(GreensFunc& green, const Walker& walker, TimeIndex time_index, Spin spin) = 0;
    };
}

#endif // MODEL_BASE_H