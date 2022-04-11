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

    class ModelBase {
        protected:
            typedef int SpaceIndex;
            typedef int TimeIndex;
            typedef int Spin;
            typedef QuantumMonteCarlo::DqmcWalker Walker;
            typedef Lattice::LatticeBase Lattice;
            typedef Eigen::MatrixXd GreensFunc;

            ModelBase() = default;

            virtual void initial(const Lattice& lattice);

            virtual double get_update_radio(SpaceIndex space_index, TimeIndex time_index, const Walker& walker);

            virtual void local_update(SpaceIndex space_index, TimeIndex time_index, const Walker& walker);

            virtual void mult_B_from_left(GreensFunc& green, TimeIndex time_index, Spin spin, const Walker& walker);
            virtual void mult_B_from_right(GreensFunc& green, TimeIndex time_index, Spin spin, const Walker& walker);
            virtual void mult_invB_from_left(GreensFunc& green, TimeIndex time_index, Spin spin, const Walker& walker);
            virtual void mult_invB_from_right(GreensFunc& green, TimeIndex time_index, Spin spin, const Walker& walker);
            virtual void mult_transB_from_left(GreensFunc& green, TimeIndex time_index, Spin spin, const Walker& walker);
    };
}

#endif // MODEL_BASE_H