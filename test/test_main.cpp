#include <iostream>

#include "lattice/lattice_base.h"
#include "lattice/square2d.h"

#include "model/model_base.h"
#include "model/repulsive_hubbard.h"

#include "measure/observable.h"
#include "measure/observable_handler.h"
#include "measure/measure_handler.h"

#include "dqmc_walker.h"
#include "dqmc.h"
#include "dqmc_initializer.h"

#include "svd_stack.h"
#include "fft_solver.h"
#include "utils/linear_algebra.hpp"
#include "utils/numerical_stable.hpp"
#include "utils/progress_bar.hpp"
#include "random.h"

#include "checkerboard/checkerboard_base.h"
#include "checkerboard/square2d.h"



// #include "random.h"
// #include "hubbard.h"

// #define EIGEN_USE_MKL_ALL
// #define EIGEN_VECTORIZE_SSE4_2
// #include <Eigen/Core>


int main() {

    // test Repulsive Hubbard

    // some params
    int ll = 2;
    double beta = 4.0;
    int lt = 80;

    int nwrap = 10;

    double hopping_t = 1.0;
    double onsite_u  = 4.0;
    double chemical_potential = 0.0;

    // fixed random seed for debug
    Utils::Random::set_seed_fixed(12345);

    Model::ModelBase* model = new Model::RepulsiveHubbard();
    Lattice::LatticeBase* lattice = new Lattice::Square2d();
    QuantumMonteCarlo::DqmcWalker* walker = new QuantumMonteCarlo::DqmcWalker();
    Measure::MeasureHandler* meas_handler = new Measure::MeasureHandler();

    // set up params
    lattice->set_space_size(ll);
    walker->set_physical_params(beta, lt);
    walker->set_stabilization_pace(nwrap);
    model->set_model_params(hopping_t, onsite_u, chemical_potential);

    // initialize modules
    QuantumMonteCarlo::DqmcInitializer::initial_modules(*lattice, *model, *walker, *meas_handler);

    model->set_bosonic_fields_to_random();

    QuantumMonteCarlo::DqmcInitializer::initial_dqmc(*lattice, *model, *walker, *meas_handler);

    // std::cout << walker->GreenttUp() << std::endl;
    // std::cout << lattice->HoppingMatrix() << std::endl;


    // test checker board
    // even lattice size and efficiency ?
    CheckerBoard::Square2d* cb = new CheckerBoard::Square2d();
    cb->set_params(ll, ll*ll, walker->TimeInterval(), hopping_t, chemical_potential);
    cb->initial();












    // // test ObservableHandler 
    // Observable::ObservableHandler* handler = new Observable::ObservableHandler();

    // std::vector<std::string> obs_list{ "filling_number", };
    // handler->initial(obs_list);

    // if ( handler->find("filling_number") ) {
    //     std::cout << "found!" << std::endl;
    //     const auto obs = handler->find_scalar("filling_number");
    //     std::cout << obs.name() << std::endl;
    // }
    // if ( handler->find("eqtime_sign") ) {
    //     std::cout << "found!" << std::endl;
    //     const auto obs = handler->find_scalar("eqtime_sign");
    //     std::cout << obs.name() << std::endl;
    // }
    // std::cout << (handler->m_eqtime_scalar_obs[0])->name() << std::endl;
    
    // Measure::MeasureHandler* meas_handler = new Measure::MeasureHandler();
    // Model::ModelBase* model = new Model::ModelBase();
    // Lattice::LatticeBase* lattice = new Lattice::Square2d();

    // (handler->m_eqtime_scalar_obs[0])->measure(*meas_handler, *model, *lattice);








    // // test utils

    // Utils::SvdStack* svd_stack = new Utils::SvdStack(4,10);
    // // Utils::FFTSolver::FFTSolver2d* solver = new Utils::FFTSolver::FFTSolver2d();

    // Eigen::MatrixXd mat = Eigen::MatrixXd::Random(4,4);
    
    // svd_stack->push(mat);
    // svd_stack->push(mat);

    // const auto& u = svd_stack->MatrixU();
    // const auto& s = svd_stack->SingularValues().asDiagonal();
    // const auto& v = svd_stack->MatrixV();
    // std::cout << (u*s*v.transpose() - mat*mat).maxCoeff() << std::endl;

    // Utils::Random::set_seed(1);









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