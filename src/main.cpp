#include "hubbard.h"
#include "detQMC.h"

#include <boost/program_options.hpp>
//#include <boost/format.hpp>

/**
 *  TODO:
 *   1. get params from command lines, using boost (done)
 *   2. equal-time measurements of momentum distribution and spin-spin correlation (done)
 *   3. bin measurements (done)
 *   4. time-displaced green function and measurements (done)
 *   5. ******** Modify command console output ******** (done)
 *   6. attractive interaction U < 0 (done)
 *   7. determine the critical temperature of superconducting transition (done)
 *   8. reweighing for doped case (done)
 *   9. read aux field configurations from input file (done)
 *   10. checkerboard decomposition (done)
 *   11. openmp parallel sampling (missing)
 *   12. new feature in standard of c++20, modify message output using std::format() (missing)
 *   13. log output (missing)
 *   14. simulate with bash script (missing)
 *   15. ...
 */


/**  The Main Program */
int main(int argc, char* argv[]) {

    /** model and controlling params */
    int ll = 4;
    int lt = 80;
    double beta = 4.0;
    double t = 1.0;
    double u = -4.0;
    double mu = 0.0;
    bool bool_checkerboard = false;

    int nwrap = 10;
    int nwarm = (int)(4 * ll * ll * beta);

    int nbin = 20;
    int nsweep = 100;
    int nBetweenBins = 10;

    std::string filename_eqtime = "../results/meas-eqtime.dat";
    std::string filename_dynamic = "../results/meas-dynamic.dat";
    bool bool_append = true;
    bool bool_display_process = true;

    bool bool_warm_up = true;
    bool bool_measure_eqtime = true;
    bool bool_measure_dynamic = true;


    /** read params from command line */
    boost::program_options::options_description opts("Program options");
    boost::program_options::variables_map vm;

    opts.add_options()
            ("help,h", "display this information")
            ("ll", boost::program_options::value<int>(&ll)->default_value(4), "spatial size of lattice, default: 4")
            ("lt", boost::program_options::value<int>(&lt)->default_value(80), "imaginary-time size of lattice, default: 80")
            ("beta", boost::program_options::value<double>(&beta)->default_value(4.0), "inverse temperature, default: 4.0")
            ("t", boost::program_options::value<double>(&t)->default_value(1.0), "hopping strength, default: 1.0")
            ("u", boost::program_options::value<double>(&u)->default_value(-4.0),
                    "interaction strength, u > 0 for repulsive and u < 0 for attractive case, default: -4.0")
            ("mu", boost::program_options::value<double>(&mu)->default_value(0.0), "chemical potential, default: 0.0")
            ("checkerboard", boost::program_options::value<bool>(&bool_checkerboard)->default_value(false), "whether to perform checkerboard break-up, default: false")
            ("nwrap", boost::program_options::value<int>(&nwrap)->default_value(10), "pace of stabilization process, default: 10")
            ("nwarm", boost::program_options::value<int>(&nwarm)->default_value((int)(4*ll*ll*beta)), "number of warmup sweeps, default: 4*ll*ll*beta")
            ("nbin", boost::program_options::value<int>(&nbin)->default_value(20), "number of bins, default: 20")
            ("nsweep", boost::program_options::value<int>(&nsweep)->default_value(100), "number of measurement sweeps in a bin, default: 100")
            ("nbetweenbins", boost::program_options::value<int>(&nBetweenBins)->default_value(10),
                    "number of sweeps between bins to avoid correlation, default: 10")
            ("app", boost::program_options::value<bool>(&bool_append)->default_value(true), "outfile mode: app or trunc, default: true")
            ("eqtime", boost::program_options::value<bool>(&bool_measure_eqtime)->default_value(true), "whether to do equal-time measurements, default: true")
            ("dynamic", boost::program_options::value<bool>(&bool_measure_dynamic)->default_value(true), "whether to do dynamic measurements, default: true")
            ("oeq", boost::program_options::value<std::string>(&filename_eqtime)->default_value("../results/meas-eqtime.dat"),
                    "output filename of equal-time data, default: ../results/meas-eqtime.dat")
            ("ody", boost::program_options::value<std::string>(&filename_dynamic)->default_value("../results/meas-dynamic.dat"),
                    "output filename of dynamic data, default: ../results/meas-dynamic.dat");

    try {
        boost::program_options::store(parse_command_line(argc, argv, opts), vm);
    }
    catch (...) {
        std::cerr << "Got undefined options from command line! "<< std::endl;
        exit(1);
    }
    boost::program_options::notify(vm);

    if (vm.count("help")) {
        std::cerr << argv[0] << std::endl;
        std::cerr << opts << std::endl;
        exit(1);
    }

    if ((!vm["ll"].defaulted() || !vm["beta"].defaulted()) && vm["nwarm"].defaulted()) {
        nwarm = 4 * vm["ll"].as<int>() * vm["ll"].as<int>() * (int)vm["beta"].as<double>();
    }


    /** DQMC Simulation */
    detQMC dqmc;

    /** Measure observable quantities over interaction strength U */

    std::vector<double> list_u = { 4.0, };

    for (auto uint : list_u) {
        bool_append = true;

        dqmc.set_model_params(ll, lt, beta, t, uint, mu, nwrap, bool_checkerboard);

        dqmc.set_Monte_Carlo_params(nwarm, nbin, nsweep, nBetweenBins);

        dqmc.set_controlling_params(bool_warm_up, bool_measure_eqtime, bool_measure_dynamic);

        dqmc.set_lattice_momentum(0.5, 0.5);

        dqmc.print_params();

        dqmc.init_measure();

        dqmc.run_QMC(bool_display_process);

        dqmc.analyse_stats();

        dqmc.print_stats();

        dqmc.file_output_tau_seq("../results/tau.dat");

        dqmc.file_output_stats_in_bins_dynamic("../results/g_bin.dat");

        dqmc.file_output_stats_dynamic("../results/dynamic.dat", false);

//        std::stringstream ss;
//        ss << std::setiosflags(std::ios::fixed) << std::setprecision(0) << ll;
//        std::string str_l = ss.str();
//        ss.str("");
//        ss << std::setiosflags(std::ios::fixed) << std::setprecision(0) << beta;
//        std::string str_beta = ss.str();
//        ss.str("");
//        filename_eqtime = "../results/eqtime_L" + str_l + "_beta" + str_beta + "_repulsive.dat";
//
//        dqmc.file_output_stats_eqtime(filename_eqtime, bool_append);
//
//        dqmc.file_output_stats_dynamic(filename_dynamic, bool_append);
    }


    /** Measure observable quantities in momentum space ( fermi surface ) */
//
//    std::vector<double> list_u = { -4.0, -2.0, -1.0, };
//
//    for (auto uint : list_u) {
//
//        dqmc.set_model_params(ll, lt, beta, t, uint, mu, nwrap, bool_checkerboard);
//        dqmc.set_Monte_Carlo_params(nwarm, nbin, nsweep, nBetweenBins);
//
//        std::stringstream ss;
//        ss << std::setiosflags(std::ios::fixed) << std::setprecision(0) << ll;
//        std::string str_l = ss.str();
//        ss.str("");
//        ss << std::setiosflags(std::ios::fixed) << std::setprecision(0) << beta;
//        std::string str_beta = ss.str();
//        ss.str("");
//        ss << std::setiosflags(std::ios::fixed) << std::setprecision(1) << uint;
//        std::string str_u = ss.str();
//        ss.str("");
//        std::string filename = "../results/fermi_surface_L" + str_l + "_beta" + str_beta + "_u" + str_u + ".dat";
//
//        for (int i = 0; i <= ll; ++i) {
//            // crystal momentum qx
//            const double qx = 1.0 - 2 * i / ll;
//
//            // time reverse symmetry: G(k, beta/2) = G(-k, beta/2)
//            for (int j = 0; j <= ll - i; ++j) {
//                // crystal momentum qy
//                const double qy = 1.0 - 2 * j / ll;
//
//                dqmc.set_lattice_momentum(qx, qy);
//
//                // warm up only one time
//                bool_warm_up = (i == 0 && j == 0);
//
//                dqmc.set_controlling_params(bool_warm_up, false, true);
//
//                dqmc.print_params();
//
//                dqmc.init_measure();
//
//                dqmc.run_QMC(bool_display_process);
//
//                dqmc.analyse_stats();
//
//                dqmc.print_stats();
//
//                bool_append = !(i == 0 && j == 0);
//                std::ofstream outfile;
//                outfile.open(filename, std::ios::out | ((bool_append)? std::ios::app : std::ios::trunc));
//                outfile << std::setiosflags(std::ios::right)
//                        << std::setw(15) << i
//                        << std::setw(15) << j
//                        << std::setw(15) << qx
//                        << std::setw(15) << qy
//                        << std::setw(15) << dqmc.dynamicMeasure.obs_mean_g_kt[ceil(lt/2)]
//                        << std::setw(15) << dqmc.dynamicMeasure.obs_err_g_kt[ceil(lt/2)]
//                        << std::endl;
//                outfile.close();
//            }
//        }
//    }

    /** Measure dynamic green's function */
//
//    dqmc.set_model_params(ll, lt, beta, t, u, mu, nwrap, bool_checkerboard);
//
//    dqmc.set_Monte_Carlo_params(nwarm, nbin, nsweep, nBetweenBins);
//
//    dqmc.set_controlling_params(bool_warm_up, false, bool_measure_dynamic);
//
//    dqmc.set_lattice_momentum(0.5, 0.5);
//
//    dqmc.print_params();
//
//    dqmc.init_measure();
//
//    dqmc.run_QMC(bool_display_process);
//
//    dqmc.analyse_stats();
//
//    dqmc.print_stats();
//
//    std::stringstream ss;
//    ss << std::setiosflags(std::ios::fixed) << std::setprecision(1) << beta;
//    std::string str_beta = ss.str();
//    ss.str("");
//    ss << std::setiosflags(std::ios::fixed) << std::setprecision(0) << ll;
//    std::string str_l = ss.str();
//    ss.str("");
//    ss << std::setiosflags(std::ios::fixed) << std::setprecision(0) << lt;
//    std::string str_lt = ss.str();
//    ss.str("");
//    ss << std::setiosflags(std::ios::fixed) << std::setprecision(1) << u;
//    std::string str_u = ss.str();
//    ss.str("");
//    std::string filename = "../results/gt_l" + str_l + "_lt" + str_lt + "_u" + str_u + "_b" + str_beta + "_k_pi2pi2.dat";
//
//    std::ofstream outfile;
//    outfile.open(filename, std::ios::out | std::ios::trunc);
//
//    outfile << std::setiosflags(std::ios::right);
//    for (int l = 1; l <= lt; ++l) {
//        outfile << std::setw(15) << l
//                << std::setw(15) << dqmc.dynamicMeasure.obs_mean_g_kt[l-1]
//                << std::setw(15) << dqmc.dynamicMeasure.obs_err_g_kt[l-1]
//                << std::endl;
//    }
//    outfile.close();


    /** Measure helicity modules over temperature T */
//
//    std::vector<double> list_beta = { 5.0, 5.0, 5.0, 8.0, };
//
//    for (auto Beta : list_beta) {
//
//        lt = (int)(Beta / 0.05);
//
//        dqmc.set_model_params(ll, lt, Beta, t, u, mu, nwrap, bool_checkerboard);
//
//        dqmc.set_Monte_Carlo_params(nwarm, nbin, nsweep, nBetweenBins);
//
//        std::stringstream ss;
//        ss << std::setiosflags(std::ios::fixed) << std::setprecision(1) << Beta;
//        std::string str_beta = ss.str();
//        ss.str("");
//        ss << std::setiosflags(std::ios::fixed) << std::setprecision(0) << ll;
//        std::string str_l = ss.str();
//        ss.str("");
//        ss << std::setiosflags(std::ios::fixed) << std::setprecision(0) << lt;
//        std::string str_lt = ss.str();
//        ss.str("");
//        ss << std::setiosflags(std::ios::fixed) << std::setprecision(1) << u;
//        std::string str_u = ss.str();
//        ss.str("");
//
//        const std::string fileConfigs = "../results/rhos_L_" + str_l + "_u_" + str_u +"/config_L_" + str_l + "_lt_" + str_lt + "_b_" + str_beta + ".dat";
//        std::ifstream infile;
//        infile.open(fileConfigs, std::ios::in);
//
//        if (!infile.is_open()) {
//            std::cerr << "fail to open file " + fileConfigs + ", start simulation with random configs." << std::endl;
//            bool_warm_up = true;
//        }
//        else {
//            infile.close();
//            dqmc.read_aux_field_configs(fileConfigs);
//            std::cerr << "old configuration is read from " + fileConfigs +", no need to warm up." << std::endl;
//            bool_warm_up = false;
//        }
//
//        dqmc.set_controlling_params(bool_warm_up, false, bool_measure_dynamic);
//
//        dqmc.set_lattice_momentum(1.0, 1.0);
//
//        dqmc.print_params();
//
//        dqmc.init_measure();
//
//        dqmc.run_QMC(bool_display_process);
//
//        dqmc.analyse_stats();
//
//        dqmc.print_stats();
//
//        dqmc.file_output_aux_field_configs(fileConfigs);
//
//        bool_append = true;
//        std::string filename = "../results/rhos_L_" + str_l + "_u_" + str_u + "/sc_rhos_L_" + str_l + "_u_" + str_u + ".dat";
//        std::string filename_bins = "../results/rhos_L_" + str_l  + "_u_" + str_u + "/bins_L_" + str_l + "_u_" + str_u + "_b_" + str_beta + ".dat";
//
//        std::ofstream outfile;
//        outfile.open(filename, std::ios::out | ((bool_append)? std::ios::app : std::ios::trunc));
//        outfile << std::setiosflags(std::ios::right)
//                << std::setw(15) << Beta
//                << std::setw(15) << 1 / Beta
//                << std::setw(15) << dqmc.dynamicMeasure.obs_mean_rho_s
//                << std::setw(15) << dqmc.dynamicMeasure.obs_err_rho_s
//                << std::endl;
//        outfile.close();
//
//        outfile.open(filename_bins, std::ios::out | std::ios::app);
//        for (int bin = 0; bin < nbin; ++bin) {
//            outfile << std::setiosflags(std::ios::right)
//                    << std::setw(15) << bin + 1
//                    << std::setw(15) << dqmc.dynamicMeasure.obs_bin_rho_s[bin]
//                    << std::endl;
//        }
//        outfile.close();
//    }


    /** Checkerboard Benchmark */

//    std::chrono::steady_clock::time_point begin_t, end_t;
//
//    Hubbard hubbard1(20, 80, 4.0, 1.0, 4.0, 0.0, 10, true);
//    Hubbard hubbard2(20, 80, 4.0, 1.0, 4.0, 0.0, 10, false);
//
//    Eigen::MatrixXd test1 = Eigen::MatrixXd::Identity(hubbard1.ls, hubbard1.ls);
//    Eigen::MatrixXd test2 = Eigen::MatrixXd::Identity(hubbard2.ls, hubbard2.ls);
//
//    begin_t = std::chrono::steady_clock::now();
//    for (int i = 0; i < 100; ++i) {
//        hubbard1.mult_B_from_left(test1, 0, +1);
//    }
//    end_t = std::chrono::steady_clock::now();
//    std::cout << (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin_t).count()/1000 << std::endl;
//
//    begin_t = std::chrono::steady_clock::now();
//    for (int i = 0; i < 100; ++i) {
//        hubbard2.mult_B_from_left(test2, 0, +1);
//    }
//    end_t = std::chrono::steady_clock::now();
//    std::cout << (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin_t).count()/1000 << std::endl;

    return 0;
}
