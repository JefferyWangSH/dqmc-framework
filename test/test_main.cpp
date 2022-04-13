#include <iostream>

#include "lattice/lattice_base.h"
#include "lattice/square2d.h"

#include "model/model_base.h"

#include "measure/observable.h"
#include "measure/observable_handler.h"

#include "dqmc_walker.h"

#include "svd_stack.h"
#include "fft_solver.h"
#include "utils/linear_algebra.hpp"
#include "utils/numerical_stable.hpp"
#include "utils/progress_bar.hpp"
#include "utils/random.hpp"



// #include "random.h"
// #include "hubbard.h"

// #define EIGEN_USE_MKL_ALL
// #define EIGEN_VECTORIZE_SSE4_2
// #include <Eigen/Core>


int main() {



    // test utils

    Utils::SvdStack* svd_stack = new Utils::SvdStack(4,10);
    // Utils::FFTSolver::FFTSolver2d* solver = new Utils::FFTSolver::FFTSolver2d();

    Eigen::MatrixXd mat = Eigen::MatrixXd::Random(4,4);
    
    svd_stack->push(mat);
    svd_stack->push(mat);

    const auto& u = svd_stack->MatrixU();
    const auto& s = svd_stack->SingularValues().asDiagonal();
    const auto& v = svd_stack->MatrixV();
    std::cout << (u*s*v.transpose() - mat*mat).maxCoeff() << std::endl;

    Utils::Random::set_seed(1);









    // // test measure

    // Observable::ObservableBase* obs = new Observable::Observable<Observable::ScalarType>();

    // Observable::Observable<Observable::ScalarType>* casted_obs = (dynamic_cast<Observable::Observable<Observable::ScalarType>*>(obs));
    // casted_obs->set_observable_name("filling");
    // std::cout << casted_obs->name() << std::endl;

    // Observable::ObservableHandler* handler = new Observable::ObservableHandler();








    // // test model
    // Model::ModelBase* model = new Model::ModelBase();

    




    // // test lattice

    // Lattice::Square2d mylattice(4);
    // mylattice.initial();

    // std::cout << mylattice.SpaceDim() << std::endl;
    // std::cout << mylattice.SpaceSize() << std::endl;
    // std::cout << mylattice.TotalSiteNum() << std::endl;

    // std::cout << mylattice.site2index({2,3}) << std::endl;
    // auto site = mylattice.index2site(11);
    // for (long unsigned int i = 0; i < site.size(); ++i) {
    //     std::cout << site[i] << std::endl;
    // }

    // std::cout << mylattice.HoppingMatrix() << std::endl;

    // Lattice::LatticeBase* lattice = new Lattice::Square2d(3);
    // lattice->initial();
    // std::cout << lattice->HoppingMatrix() << std::endl;

    return 0;
}