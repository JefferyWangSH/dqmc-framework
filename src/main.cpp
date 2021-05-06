#include "hubbard.h"
#include "detQMC.h"

#include <boost/program_options.hpp>

/**
 *  todo:
 *   1. get params from command lines, using boost (done)
 *   2. equal-time measurements of momentum distribution and spin-spin correlation (done)
 *   3. bin measurements (done)
 *   4. time-displaced green function and measurements (done)
 *   5. Stochastic Analytic Continuation (SAC) to obtain fermion spectrum function (missing)
 *   6. Check-board decomposition (missing)
 *   7. ******** Modify command console output ******** (done)
 *   8. attractive interaction U < 0 (done)
 *   9. openmp parallel programming (missing)
 *   10. determine the critical temperature of superconducting transition (missing)
 *   11. read aux field configurations from input file (missing)
 *   12. ...
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

    int nwrap = 10;
    int nwarm = (int)(4 * ll * ll * beta);

    int nbin = 20;
    int nsweep = 100;
    int nBetweenBins = 10;

    std::string filename_eqtime = "../results/meas-eqtime.txt";
    std::string filename_dynamic = "../results/meas-dynamic.txt";
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
            ("nwrap", boost::program_options::value<int>(&nwrap)->default_value(10), "pace of stabilization process, default: 10")
            ("nwarm", boost::program_options::value<int>(&nwarm)->default_value((int)(4*ll*ll*beta)), "number of warmup sweeps, default: 4*ll*ll*beta")
            ("nbin", boost::program_options::value<int>(&nbin)->default_value(20), "number of bins, default: 20")
            ("nsweep", boost::program_options::value<int>(&nsweep)->default_value(100), "number of measurement sweeps in a bin, default: 100")
            ("nbetweenbins", boost::program_options::value<int>(&nBetweenBins)->default_value(10),
                    "number of sweeps between bins to avoid correlation, default: 10")
            ("app", boost::program_options::value<bool>(&bool_append)->default_value(true), "outfile mode: app or trunc, default: true")
            ("eqtime", boost::program_options::value<bool>(&bool_measure_eqtime)->default_value(true), "whether to do equal-time measurements, default: true")
            ("dynamic", boost::program_options::value<bool>(&bool_measure_dynamic)->default_value(true), "whether to do dynamic measurements, default: true")
            ("oeq", boost::program_options::value<std::string>(&filename_eqtime)->default_value("../results/meas-eqtime.txt"),
                    "output filename of equal-time data, default: ../results/meas-eqtime.txt")
            ("ody", boost::program_options::value<std::string>(&filename_dynamic)->default_value("../results/meas-dynamic.txt"),
                    "output filename of dynamic data, default: ../results/meas-dynamic.txt");

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

    std::vector<double> list_u = { -4.0, };

    for (auto uint : list_u) {
        bool_append = false;

        dqmc.set_Model_Params(ll, lt, beta, t, uint, mu, nwrap);

        dqmc.set_MC_Params(nwarm, nbin, nsweep, nBetweenBins);

        dqmc.set_bool_Params(bool_warm_up, bool_measure_eqtime, bool_measure_dynamic);

        dqmc.set_Momentum_q(M_PI / 2, M_PI / 2);

        dqmc.printParams();

        dqmc.initialMeasure();

        dqmc.runQMC(bool_display_process);

        dqmc.analyseStats();

        dqmc.printStats();

        dqmc.output_Stats_eqtime(filename_eqtime, bool_append);

        dqmc.output_Stats_dynamic(filename_dynamic, bool_append);
    }


    /** Measure observable quantities in momentum space */
    /*
    dqmc.set_MC_Params(nwarm, nbin, nsweep, nBetweenBins);
    dqmc.set_Model_Params(ll, lt, beta, t, u, mu, nwrap);

    std::stringstream ss;
    ss << std::setiosflags(std::ios::fixed) << std::setprecision(0) << ll;
    std::string str_l = ss.str();
    ss.str("");
    ss << std::setiosflags(std::ios::fixed) << std::setprecision(0) << beta;
    std::string str_beta = ss.str();
    ss.str("");
    ss << std::setiosflags(std::ios::fixed) << std::setprecision(1) << u;
    std::string str_u = ss.str();
    ss.str("");
    std::string filename = "../results/fermi_surface_L" + str_l + "_beta" + str_beta + "_u" + str_u + ".txt";

    for (int i = 0; i <= ll; ++i) {
        // crystal momentum qx
        const double qx = M_PI - 2 * M_PI / ll * i;

        // time reverse symmetry: G(k, beta/2) = G(-k, beta/2)
        for (int j = 0; j <= ll - i; ++j) {
            // crystal momentum qy
            const double qy = M_PI - 2 * M_PI / ll * j;

            dqmc.set_Momentum_q(qx, qy);

            // warm up only one time
            bool_warm_up = (i == 0 && j == 0);

            dqmc.set_bool_Params(bool_warm_up, false, true);

            dqmc.printParams();

            dqmc.initialMeasure();

            dqmc.runQMC(bool_display_process);

            dqmc.analyseStats();

            dqmc.printStats();

            bool_append = !(i == 0 && j == 0);
            std::ofstream outfile;
            outfile.open(filename, std::ios::out | ((bool_append)? std::ios::app : std::ios::trunc));
            outfile << std::setiosflags(std::ios::right)
                    << std::setw(15) << i
                    << std::setw(15) << j
                    << std::setw(15) << qx
                    << std::setw(15) << qy
                    << std::setw(15) << dqmc.dynamicMeasure.obs_mean_g_kt[ceil(lt/2)]
                    << std::setw(15) << dqmc.dynamicMeasure.obs_err_g_kt[ceil(lt/2)]
                    << std::endl;
            outfile.close();
        }
    }
    */


    /** Measure observable quantities over temperature T */
    /*
    for (double T = 0.2; T <= 0.2; T += 0.2) {
        bool_append = true;

        dqmc.set_Model_Params(ll, lt, 1/T, t, u, mu, nwrap);

        dqmc.set_MC_Params(nwarm, nbin, nsweep, nBetweenBins);

        dqmc.set_bool_Params(bool_warm_up, bool_measure_eqtime, bool_measure_dynamic);

        dqmc.set_Momentum_q(M_PI, M_PI);

        dqmc.printParams();

        dqmc.initialMeasure();

        dqmc.runQMC(false);

        dqmc.analyseStats();

        dqmc.printStats();

        dqmc.output_Stats_eqtime(filename_eqtime, bool_append);

        dqmc.output_Stats_dynamic(filename_dynamic, bool_append);
    }
    */

    return 0;
}
