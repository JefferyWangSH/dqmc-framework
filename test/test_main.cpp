#include <iostream>

#include "lattice/lattice_base.h"
#include "lattice/square2d.h"

#include "model/model_base.h"
#include "model/repulsive_hubbard.h"
#include "model/repulsive_hubbard_cbsquare2d.h"

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

#include <chrono>

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

    Model::ModelBase* model = new Model::RepulsiveHubbardCbSquare2d();
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

    walker->sweep_from_0_to_beta(*model);
    walker->sweep_from_beta_to_0(*model);
    walker->sweep_for_dynamic_greens(*model);
    walker->sweep_from_beta_to_0(*model);

    std::cout << walker->GreenttUp() << std::endl;
    std::cout << std::endl;
    std::cout << walker->Greent0Up() << std::endl;


    // todo: test checkerboard
    // even lattice size (ok) and efficiency (ok)
    
    // todo: trans mult ( finished, trans is actually unnecessary, just use expK )
    // todo: Vmat mult benchmark ( two-times faster! )
    // !!!



    // Model::ModelBase* model_cb = new Model::RepulsiveHubbardCbSquare2d();
    // Model::ModelBase* model_direct = new Model::RepulsiveHubbard();
    // model_cb->set_model_params(hopping_t, onsite_u, chemical_potential);
    // model_direct->set_model_params(hopping_t, onsite_u, chemical_potential);

    // Lattice::LatticeBase* lattice = new Lattice::Square2d();
    // QuantumMonteCarlo::DqmcWalker* walker = new QuantumMonteCarlo::DqmcWalker();
    // Measure::MeasureHandler* meas_handler = new Measure::MeasureHandler();
    // lattice->set_space_size(ll);
    // walker->set_physical_params(beta, lt);
    // walker->set_stabilization_pace(nwrap);
    
    // QuantumMonteCarlo::DqmcInitializer::initial_modules(*lattice, *model_cb, *walker, *meas_handler);
    // QuantumMonteCarlo::DqmcInitializer::initial_modules(*lattice, *model_direct, *walker, *meas_handler);

    // Utils::Random::set_seed_fixed(12345);
    // model_cb->set_bosonic_fields_to_random();
    // Utils::Random::set_seed_fixed(12345);
    // model_direct->set_bosonic_fields_to_random();

    // Eigen::MatrixXd mat_cb = Eigen::MatrixXd::Identity(ll*ll, ll*ll);
    // Eigen::MatrixXd mat_direct = mat_cb;

    // const int num_mult = 5000;

    // std::chrono::steady_clock::time_point begin_t{}, end_t{};

    // begin_t = std::chrono::steady_clock::now();
    // for (auto i = 0; i < num_mult; ++i) {
    //     model_direct->mult_B_from_left(mat_direct, 0, 1);
    // }
    // end_t = std::chrono::steady_clock::now();
    // std::cout << "direct : " << std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin_t).count() << std::endl;
    
    // begin_t = std::chrono::steady_clock::now();
    // for (auto i = 0; i < num_mult; ++i) {
    //     model_cb->mult_B_from_left(mat_cb, 0, 1);
    // }
    // end_t = std::chrono::steady_clock::now();
    // std::cout << "cb : " << std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin_t).count() << std::endl;

    // // std::cout << (mat_direct - mat_cb).maxCoeff() << std::endl;









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