#ifndef DQMC_HUBBARD_MEASURE_H
#define DQMC_HUBBARD_MEASURE_H
#pragma once

/**
  *  This head file includes module for both equal-time and dynamical measuring.
  *  Class: Measure::Measure
  *  Measuring:
  *   1. Filling number <n>
  *   2. Double occupancy D = < n_up*n_dn >
  *   3. Single particle kinetic energy
  *   4. Distributions of electrons in momentum space
  *   5. Local spin correlation, magnetization C(0,0) = < (n_up - n_dn)^2 >
  *   6. Spin density structure factor (SDW)
  *   7. Charge density structure factor (CDW)
  * 
  *   8. Dynamical green's function of imaginary time: G(k, tau) = < c(k, tau) * c^+(k, 0) >
  *   9. Superfluid stiffness \rho_s of superconducting: \rho_s = (\Gamma_L - \Gamma_T) / 4
  *  10. Density of states in imaginary time space: N(tau) = 1/N * \sum_{i} G(tau, 0)_{ii}
  */

#include <memory>
#include <vector>
#include "observable.h"
#include "measure_container.h"


// forward declaration
namespace Model { class Hubbard; }
namespace Simulation { class DetQMC; }

namespace Measure {
    
    class Measure {
    private:
        // number of bins
        int _nbin{};

        // list of observables to be measured
        std::vector<std::string> _obs_list{};

        // container for physical observables
        Container _container;

        // lattice momentum q
        Eigen::VectorXd q = Eigen::VectorXd::Zero(2);

        friend class Methods;
        friend class GatherMPI;
        friend class Simulation::DetQMC;

    public:

        /* (de)construct functions */
        Measure()  = default;
        ~Measure() = default;

        /* interfaces */
        int nbin() const;
        bool is_eqtime_measure() const;
        bool is_dynamic_measure() const;
        std::vector<Observable<double>> obs_list_double() const;
        std::vector<Observable<Eigen::VectorXd>> obs_list_vector() const;
        std::vector<Observable<Eigen::MatrixXd>> obs_list_matrix() const;
        Observable<double> sign_eqtime() const;
        Observable<double> sign_dynamic() const;

        bool find(std::string obs_name) const; 
        Observable<double> find_double_obs(std::string obs_name) const;
        Observable<Eigen::VectorXd> find_vector_obs(std::string obs_name) const;
        Observable<Eigen::MatrixXd> find_matrix_obs(std::string obs_name) const;

        /* set up the size of bins */
        void set_size_of_bin(const int &nbin);

        /* set up list of observables to be measured */
        void set_observable_list(const std::vector<std::string> &obs_list);

        /* set up lattice momentum q */
        void set_lattice_momentum(const Eigen::VectorXd &q);

        /* initialize and prepare for measurements */
        void initial(const Model::Hubbard &hubbard);

        /* clear temporary data */
        void clear_temporary();

        /* do the measurements */
        void eqtime_measure(const Model::Hubbard &hubbard);
        void dynamic_measure(const Model::Hubbard &hubbard);

        /* normalize data from scratch */
        void normalize_stats();

        /* bin collection of data */
        void write_stats_to_bins(int bin);

        /* analyse statistics by calculating means and errors */
        void analyse_stats();

    };

} // namespace Measure

#endif //DQMC_HUBBARD_MEASURE_H
