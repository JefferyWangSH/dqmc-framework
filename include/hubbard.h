#ifndef DQMC_HUBBARD_HUBBARD_H
#define DQMC_HUBBARD_HUBBARD_H
#pragma once

/**
  *  This head file includes hubbard class
  *  which is defined for the DQMC simulation of 2d fermion Hubbard model
  *  Model setup and pivotal Monte Carlo updating algorithms involved.
  */

#include <memory>
#include <vector>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include "checker_board.h"
#include "svd_stack.h"

namespace Measure { class Measure; class Methods; }
namespace Simulation { class DetQMC; }
namespace FileOutput { 
    void file_output_tau(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode); 
    void file_output_aux_field(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode); 
}
namespace ScreenOutput {
    void screen_output_params(const int &world_size, const Simulation::DetQMC &dqmc); 
    void screen_output_end_info(const Simulation::DetQMC &dqmc);
}


namespace Model {

    class Hubbard {
    private:
        // model params
        int ll{4}, ls{16}, lt{80};
        double beta{4.0}, dtau{0.1};
        double t{1.0}, u_int{4.0}, mu{0.0}, alpha{0.0};

        bool is_attractive_u{true};
        bool is_eqtime_measure{true};
        bool is_dynamic_measure{true};

        int nwrap{10};
        int current_tau{0};

        double config_sign{0.0};
        std::unique_ptr<std::vector<double>> vec_config_sign{};

        double max_wrap_error_equal{0.0};
        double max_wrap_error_dynamic{0.0};

        // aux boson fields s
        std::unique_ptr<Eigen::MatrixXd> s{};

        // checkerboard class for hopping matrix fabrication
        bool is_checkerboard{true};
        std::unique_ptr<CheckerBoard::CheckerBoard> checkerboard{};

        // equal-time greens function for both spin up and down states
        // critical quantities in DQMC simulation
        std::unique_ptr<Eigen::MatrixXd> green_tt_up{}, green_tt_dn{};
        std::unique_ptr<std::vector<Eigen::MatrixXd>> vec_green_tt_up{}, vec_green_tt_dn{};

        // time-displaced greens function for dynamic measurements
        // Matsubara greens function: G(\tau, 0) and G(0, \tau).
        // Gij(\tau, 0) = < ci(\tau) * cj^+ (0) >
        std::unique_ptr<Eigen::MatrixXd> green_t0_up{}, green_t0_dn{};
        std::unique_ptr<std::vector<Eigen::MatrixXd>> vec_green_t0_up{}, vec_green_t0_dn{};

        // Gij(0, \tau) = - < cj^+(\tau) * ci(0) >
        std::unique_ptr<Eigen::MatrixXd> green_0t_up{}, green_0t_dn{};
        std::unique_ptr<std::vector<Eigen::MatrixXd>> vec_green_0t_up{}, vec_green_0t_dn{};

        // auxiliary SvdStack class for numerical stabilization
        std::unique_ptr<SvdStack> stack_left_up{};
        std::unique_ptr<SvdStack> stack_left_dn{};
        std::unique_ptr<SvdStack> stack_right_up{};
        std::unique_ptr<SvdStack> stack_right_dn{};

        // using unique_ptr may sacrifice performance while memory-saving and safe
        // compared with creating variables directly

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Hubbard() = default;
        ~Hubbard() = default;

        /** set up params */
        void set_model_params(int ll, int lt, double beta, double t, double u_int, double mu);
        void set_bool_params(bool is_eqtime_measure, bool is_dynamic_measure, bool is_checkerboard);
        void set_stabilization_pace(int nwrap);

        /** initialization */
        void initial();

        /** sweep the space-time lattice from 0 to beta */
        void sweep_0_to_beta();

        /** sweep the space-time lattice from beta to 0 */
        void sweep_beta_to_0();

        /** sweep from beta to 0 to calculate time-displaced (dynamical) green functions */
        void sweep_0_to_beta_dynamic();

        // friend classes and functions
        friend class Simulation::DetQMC;
        friend class CheckerBoard::CheckerBoard;
        friend class Measure::Measure;
        friend class Measure::Methods;
        friend void FileOutput::file_output_tau(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode);
        friend void FileOutput::file_output_aux_field(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode); 
        friend void ScreenOutput::screen_output_params(const int &world_size, const Simulation::DetQMC &dqmc);
        friend void ScreenOutput::screen_output_end_info(const Simulation::DetQMC &dqmc);

    private:
        /** (de)allocate memory for params */
        void allocate();
        void deallocate();

        /** randomly initialize aux field */
        void init_field_to_random();

        /** initialize udv stacks for sweep use */
        void init_stacks();

//        /** compute B matrix with slice l and spin sigma given */
//        Eigen::MatrixXd make_B_l(int l, int sigma);

        /** multiply B matrix in place */
        void mult_B_from_left(Eigen::MatrixXd &A, int l, int sigma);

        void mult_B_from_right(Eigen::MatrixXd &A, int l, int sigma);

        void mult_invB_from_left(Eigen::MatrixXd &A, int l, int sigma);

        void mult_invB_from_right(Eigen::MatrixXd &A, int l, int sigma);

        void mult_transB_from_left(Eigen::MatrixXd &A, int l, int sigma);

        /** update the aux field at time slice l with Metropolis algorithm */
        void metropolis_update(int l);

        /** propagate the green's function from l to l+1 */
        void wrap_0_to_beta(int l);

        /** propagate the green's function from l to l-1 */
        void wrap_beta_to_0(int l);
    };

}


#endif //DQMC_HUBBARD_HUBBARD_H
