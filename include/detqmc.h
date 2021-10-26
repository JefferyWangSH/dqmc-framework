#ifndef DQMC_HUBBARD_DETQMC_H
#define DQMC_HUBBARD_DETQMC_H
#pragma once

/**
  *  This head file includes detQMC class
  *  for the Monte Carlo simulation of half-filled Hubbard model.
  *  including: MC sampling (measuring) and output of statistical results.
  */

class SvdStack;
namespace Model { class Hubbard; }
namespace Measure { class EqtimeMeasure; }
namespace Measure { class DynamicMeasure; }
namespace Simulation { class DetQMC; }
namespace ScreenOutput { 
    void screen_output_params(const int &world_size, const Simulation::DetQMC &dqmc); 
    void screen_output_end_info(const Simulation::DetQMC &dqmc);
}

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include <chrono>


namespace Simulation {

    class DetQMC {
    private:
        // model parameters
        Model::Hubbard *hubb{};

        int nwrap{10}, nwarm{300};
        int nbin{20}, nsweep{100}, n_between_bins{10};

        // bool parameters to control thermalization and measurements
        bool bool_warm_up = true;
        bool bool_measure_eqtime = true;
        bool bool_measure_dynamic = true;

        // lattice momentum q
        Eigen::VectorXd q = Eigen::VectorXd::Zero(2);

        // time cost of one single measuring process
        std::chrono::steady_clock::time_point begin_t{}, end_t{};

    public:
        // for equal-time measurements
        Measure::EqtimeMeasure *EqtimeMeasure{};

        // for time-displaced (dynamical) measurements
        Measure::DynamicMeasure *DynamicMeasure{};

        // friend function
        friend void ScreenOutput::screen_output_params(const int &world_size, const Simulation::DetQMC &dqmc);
        friend void ScreenOutput::screen_output_end_info(const Simulation::DetQMC &dqmc);

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        DetQMC() = default;

        ~DetQMC();

        /** functions to set up params and initialization before DQMC calculation starts */
        /* set up model parameters */
        void set_model_params(int ll, int lt, double beta, double t, double Uint, double mu, int nwrap, bool is_checkerboard);

        /* set up parameters for Monte Carlo simulation */
        void set_Monte_Carlo_params(int nwarm, int nbin, int nsweep, int n_between_bins);

        /* set up bool parameters */
        void set_controlling_params(bool bool_warm_up, bool bool_measure_eqtime, bool bool_measure_dynamic);

        /* set up lattice momentum q for momentum measurements */
        void set_lattice_momentum(double qx, double qy);

        /* read aux field configurations from input file */
        void read_aux_field_configs(const std::string &filename);


        /** Critical functions for DQMC calculations, including Monte Carlo updates and measurements **/
        /* initialization, especially preparing for measuring */
        void init_measure();

        /* run a dqmc simulation */
        void run(bool bool_display_process);

        /* analyse statistics from simulation */
        void analyse_stats() const;

    private:
        /* process of Monte Carlo sweeping, and do the measurements if needed */
        void sweep_back_and_forth(bool bool_eqtime, bool bool_dynamic) const;
    };

}


#endif //DQMC_HUBBARD_DETQMC_H
