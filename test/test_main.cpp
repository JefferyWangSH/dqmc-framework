#include <iostream>

#include "lattice/lattice_base.h"
#include "lattice/square.h"

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
#include "random.h"

#include "checkerboard/checkerboard_base.h"
#include "checkerboard/square.h"

#include <chrono>

#include <unistd.h>
#include "utils/progressbar.hpp"

// #include "random.h"
// #include "hubbard.h"

// #define EIGEN_USE_MKL_ALL
// #define EIGEN_VECTORIZE_SSE4_2
// #include <Eigen/Core>


int main() {

    // test Repulsive Hubbard

    // some params
    int ll = 4;
    double beta = 8.0;
    int lt = 160;

    int nwrap = 10;

    double hopping_t = 1.0;
    double onsite_u  = 4.0;
    double chemical_potential = 0.0;

    int sweeps_warmup = 512;
    int bin_num = 20;
    int bin_size = 100;
    int sweeps_between_bins = 20;

    std::vector<std::string> obs_list = { 
                                          "filling_number", 
                                          "double_occupancy",
                                          "kinetic_energy",
                                          "local_spin_corr",
                                          "greens_functions", 
                                          };

    // fixed random seed for debug
    Utils::Random::set_seed_fixed(12345);

    Model::ModelBase* model = new Model::RepulsiveHubbard();
    Lattice::LatticeBase* lattice = new Lattice::Square();
    QuantumMonteCarlo::DqmcWalker* walker = new QuantumMonteCarlo::DqmcWalker();
    Measure::MeasureHandler* meas_handler = new Measure::MeasureHandler();

    // set up params
    lattice->set_lattice_params({ll,ll});
    // lattice should be initialized once lattice params have been set.
    lattice->initial();
    walker->set_physical_params(beta, lt);
    walker->set_stabilization_pace(nwrap);
    model->set_model_params(hopping_t, onsite_u, chemical_potential);
    meas_handler->set_measure_params(sweeps_warmup, bin_num, bin_size, sweeps_between_bins);
    meas_handler->set_observables(obs_list);

    // make sure that the lattice module has been initialized
    if ( lattice->InitialStatus() ) {
        QuantumMonteCarlo::DqmcInitializer::set_measured_momentum(*meas_handler, lattice->GammaPointIndex());
        QuantumMonteCarlo::DqmcInitializer::set_measured_momentum_list(*meas_handler, lattice->kStarsIndex());
    }

    // initialize modules
    // lattice module has been initialized before, so in this function 
    // it is provided for the initialization of other modules
    QuantumMonteCarlo::DqmcInitializer::initial_modules(*lattice, *model, *walker, *meas_handler);

    // using checkerboard break-up
    // CheckerBoard::CheckerBoardBase* checkerboard = new CheckerBoard::Square();
    // QuantumMonteCarlo::DqmcInitializer::initial_modules(*lattice, *model, *walker, *meas_handler, *checkerboard);

    model->set_bosonic_fields_to_random();

    QuantumMonteCarlo::DqmcInitializer::initial_dqmc(*lattice, *model, *walker, *meas_handler);

    // walker->sweep_from_0_to_beta(*model);
    // walker->sweep_from_beta_to_0(*model);
    // walker->sweep_for_dynamic_greens(*model);
    // walker->sweep_from_beta_to_0(*model);

    // QuantumMonteCarlo::Dqmc::sweep_forth_and_back(*walker, *model, *lattice, *meas_handler);

    // std::chrono::steady_clock::time_point begin_t{}, end_t{};
    // begin_t = std::chrono::steady_clock::now();


    QuantumMonteCarlo::Dqmc::show_progress_bar( true );
    QuantumMonteCarlo::Dqmc::progress_bar_format( 70, '=', ' ' );

    QuantumMonteCarlo::Dqmc::thermalize(*walker, *model, *lattice, *meas_handler);
    // std::cout << QuantumMonteCarlo::Dqmc::timer() << std::endl;

    // int loop = 1e3;
    // for (int i = 0; i < loop; ++i) {
    //     walker->wrap_from_0_to_beta(*model, 0);
    //     walker->wrap_from_beta_to_0(*model, 1);
    // }


    // end_t = std::chrono::steady_clock::now();
    // std::cout << "warm-up : " << std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin_t).count() << std::endl;


    QuantumMonteCarlo::Dqmc::measure(*walker, *model, *lattice, *meas_handler);
    // std::cout << QuantumMonteCarlo::Dqmc::timer() << std::endl;
    QuantumMonteCarlo::Dqmc::analyse(*meas_handler);

    // if (meas_handler->find("filling_number")) {
    //     auto obs = meas_handler->find_scalar("filling_number");
    //     std::cout << obs.name() << "  " << obs.mean_value() << "  " << obs.error_bar() << std::endl;
    // }

    // std::cout << meas_handler->BinsNum() << std::endl;
    // std::cout << meas_handler->BinsSize() << std::endl;

    std::cout << " wrap error :  " << walker->WrapError() << std::endl;

    if (meas_handler->find("equaltime_sign")) {
        const auto obs = meas_handler->find<Observable::ScalarObs>("equaltime_sign");
        std::cout << obs.name() << "  " << obs.mean_value() << "  " << obs.error_bar() << std::endl;
    }

    if (meas_handler->find("dynamic_sign")) {
        auto obs = meas_handler->find<Observable::ScalarObs>("dynamic_sign");
        std::cout << obs.name() << "  " << obs.mean_value() << "  " << obs.error_bar() << std::endl;
    }
      
    if (meas_handler->find("filling_number")) {
        auto obs = meas_handler->find<Observable::ScalarObs>("filling_number");
        std::cout << obs.name() << "  " << obs.mean_value() << "  " << obs.error_bar() << std::endl;
    }

    if (meas_handler->find("double_occupancy")) {
        auto obs = meas_handler->find<Observable::ScalarObs>("double_occupancy");
        std::cout << obs.name() << "  " << obs.mean_value() << "  " << obs.error_bar() << std::endl;
    }

    if (meas_handler->find("kinetic_energy")) {
        auto obs = meas_handler->find<Observable::ScalarObs>("kinetic_energy");
        std::cout << obs.name() << "  " << obs.mean_value() << "  " << obs.error_bar() << std::endl;
    }

    if (meas_handler->find("local_spin_corr")) {
        auto obs = meas_handler->find<Observable::ScalarObs>("local_spin_corr");
        std::cout << obs.name() << "  " << obs.mean_value() << "  " << obs.error_bar() << std::endl;
    }

    if (meas_handler->find("greens_functions")) {
        auto obs = meas_handler->find<Observable::MatrixObs>("greens_functions");
        std::cout << obs.name() << std::endl;
        for (int t = 0; t < walker->TimeSize(); ++t) {
            std::cout << t << "     " 
                      << obs.mean_value()(1,t) << "      " 
                      << obs.error_bar()(1,t) 
                      << std::endl;
        }
    }

    // std::cout << lattice->Index2Momentum(5) << std::endl;

    // std::cout << walker->GreenttUp() << std::endl;
    // std::cout << std::endl;
    // std::cout << walker->Greent0Up() << std::endl;
    // std::cout << std::endl;
    // std::cout << lattice->HoppingMatrix() << std::endl;

    // std::cout << std::endl;
    // std::cout << meas_handler->isWarmUp() << std::endl;
    // std::cout << meas_handler->isEqualTime() << std::endl;
    // std::cout << meas_handler->isDynamic() << std::endl;
    // std::cout << meas_handler->find_scalar("filling_number").name() << std::endl;
    // std::cout << meas_handler->find_matrix("greens_functions").name() << std::endl;











    // // test lattice momentum
    // int ll = 4;

    // Lattice::LatticeBase* lattice = new Lattice::Square();

    // lattice->set_lattice_params({ll,ll});
    // lattice->initial();

    // std::cout << lattice->Displacement(3,2) << std::endl; 
    // std::cout << lattice->Displacement(2,4) << std::endl; 

    // // std::cout << lattice->Index2Site(1) << std::endl;
    // // std::cout << lattice->Index2Momentum(lattice->GammaPointIndex()) << std::endl;
    // // std::cout << lattice->Index2Momentum(lattice->XPointIndex()) << std::endl;
    // // std::cout << lattice->Index2Momentum(lattice->MPointIndex()) << std::endl;
    
    // // std::cout << lattice->m_index2site_table << std::endl;
    // // std::cout << lattice->m_index2momentum_table << std::endl;

    // // std::cout << lattice->m_fourier_factor_table << std::endl;
    // // std::cout << lattice->FourierFactor(1,1) << std::endl;

    // // for (auto i : lattice->m_k_stars_index) {
    // //     std::cout << i << std::endl;
    // // }

    // // std::cout << lattice->m_gamma_point_index << std::endl;
    // // std::cout << lattice->m_x_point_index << std::endl;
    // // std::cout << lattice->m_m_point_index << std::endl;

    // // for (auto i : lattice->m_delta_line_index ) {
    // //     std::cout << i << std::endl;
    // // }
    // // std::cout << std::endl;

    // // for (auto i : lattice->m_z_line_index ) {
    // //     std::cout << i << std::endl;
    // // }
    // // std::cout << std::endl;

    // // for (auto i : lattice->m_sigma_line_index ) {
    // //     std::cout << i << std::endl;
    // // }
    // // std::cout << std::endl;
    
    // // for (auto i : lattice->m_gamma2x2m2gamma_loop_index ) {
    // //     std::cout << i << std::endl;
    // // }
    // // std::cout << std::endl;
    
















    // // test progress bar
    // const int total = 10000;

    // /*
    //  * Define a progress bar that has a total of 10000,
    //  * a width of 70, shows `#` to indicate completion
    //  * and a dash '-' for incomplete
    //  */
    // progresscpp::ProgressBar progressBar(total, 70, '#', '-');

    // for (int i = 0; i < total; i++) {
        

    //     usleep(200); // simulate work

    //     ++progressBar; // record the tick

    //     // display the bar only at certain steps
    //     // if (i % 10 == 0)
        
    //     std::cout << "hello "; progressBar.display();
    // }

    // // tell the bar to finish
    // std::cout << "hello "; progressBar.done();

    // std::cout << "Done!" << std::endl;




    








    // test lattice


















    // todo: test checkerboard
    // even lattice size (ok) and efficiency (ok)
    
    // todo: trans mult ( finished, trans is actually unnecessary, just use expK )
    // todo: Vmat mult benchmark ( two-times faster! )
    // !!!


    // todo: functional ptr to mult_expK in Model (ok!)
    // ( combine std::function with std::bind to wrap the member function )
    // !!




    // Model::ModelBase* model_cb = new Model::RepulsiveHubbard();
    // Model::ModelBase* model_direct = new Model::RepulsiveHubbard();
    // model_cb->set_model_params(hopping_t, onsite_u, chemical_potential);
    // model_direct->set_model_params(hopping_t, onsite_u, chemical_potential);

    // Lattice::LatticeBase* lattice = new Lattice::Square();
    // QuantumMonteCarlo::DqmcWalker* walker = new QuantumMonteCarlo::DqmcWalker();
    // Measure::MeasureHandler* meas_handler = new Measure::MeasureHandler();
    // CheckerBoard::CheckerBoardBase* checkerboard = new CheckerBoard::Square();

    // lattice->set_lattice_params({ll,ll});
    // walker->set_physical_params(beta, lt);
    // walker->set_stabilization_pace(nwrap);
    // meas_handler->set_measure_params(sweeps_warmup, bin_num, bin_size, sweeps_between_bins);
    // meas_handler->set_observables(obs_list);
    
    // QuantumMonteCarlo::DqmcInitializer::initial_modules(*lattice, *model_cb, *walker, *meas_handler, *checkerboard);
    // QuantumMonteCarlo::DqmcInitializer::initial_modules(*lattice, *model_direct, *walker, *meas_handler);

    // Utils::Random::set_seed_fixed(12345);
    // model_cb->set_bosonic_fields_to_random();
    // Utils::Random::set_seed_fixed(12345);
    // model_direct->set_bosonic_fields_to_random();

    // Eigen::MatrixXd mat_cb = Eigen::MatrixXd::Identity(ll*ll, ll*ll);
    // Eigen::MatrixXd mat_direct = mat_cb;

    // const int num_mult = 1e4;

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

    // std::cout << (mat_direct - mat_cb).maxCoeff() << std::endl;









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