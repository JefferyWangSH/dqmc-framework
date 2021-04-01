#ifndef DQMC_HUBBARD_OPTIONS_H
#define DQMC_HUBBARD_OPTIONS_H

/**
 *  This head file includes subroutines
 *  to help read simulating params from command line.
 */

#include <map>
#include <algorithm>

int option2int(const std::string& option) {
    if ( option == "ll") { return 0;}
    if ( option == "lt") { return 1;}
    if ( option == "beta") { return 2;}
    if ( option == "t") { return 3;}
    if ( option == "u") { return 4;}
    if ( option == "mu") { return 5;}
    if ( option == "nwrap") { return 6;}
    if ( option == "nwarm") { return 7;}
    if ( option == "nbin") { return 8;}
    if ( option == "nsweep") { return 9;}
    if ( option == "nBetweenBins") { return 10;}
    if ( option == "app") { return 11;}
    if ( option == "eqtime") { return 12;}
    if ( option == "dynamic") { return 13;}
    if ( option == "oeq") { return 14;}
    if ( option == "ody") { return 15;}
    if ( option == "h") { return 16;}
    if ( option == "?") { return 17;}
    return -1;
}

double str2double(const std::string& str) {
    if (!str.empty()) {
        double output;
        std::stringstream sstream;
        sstream << str;
        sstream >> output;
        return output;
    }
    else { return 0;}
}

int str2int(const std::string& str) {
    if (!str.empty()) {
        int output;
        std::stringstream sstream;
        sstream << str;
        sstream >> output;
        return output;
    }
    else { return 0;}
}

// the usage
void usage(char *program) {
    std::cerr << "\nProgram Usage:\n"
              << program << std::endl
              << "   [ -ll  4 ]                             (spatial) size of lattice, default: 4\n"
              << "   [ -lt  80 ]                            (imaginary time) size of lattice, default: 80\n"
              << "   [ -beta  4.0 ]                         inverse temperature, default: 4.0\n"
              << "   [ -t  1.0 ]                            hopping strength, default: 1.0\n"
              << "   [ -u  4.0 ]                            interaction strength, default: 4.0\n"
              << "   [ -mu  0.0 ]                           chemical potential, default: 0.0\n"
              << "   [ -nwrap  10 ]                         pace of stabilization process, default: 10\n"
              << "   [ -nwarm  256 ]                        number of warmup sweeps, default: 256\n"
              << "   [ -nbin  20 ]                          number of bins, default: 20\n"
              << "   [ -nsweep  100 ]                       number of measurement sweeps in a bin, default: 100\n"
              << "   [ -nBetweenBins  10 ]                  number of sweeps between bins to avoid correlation, default: 10\n"
              << "   [ -app  true ]                         outfile mode: app or trunc, default: true\n"
              << "   [ -eqtime  true ]                      whether to do equal-time measurements, default: true\n"
              << "   [ -dynamic  true ]                     whether to do dynamic measurements, default: true\n"
              << "   [ -oeq  ../results/meas-eqtime.txt ]   output filename, default: ../results/meas-eqtime.txt\n"
              << "   [ -ody  ../results/meas-dynamic.txt ]  output filename, default: ../results/meas-dynamic.txt\n"
              << "   [ -h / -? ]                            display this information\n"
              << std::endl;
    exit(1);
}

// function to read and reload parameters from the command line
void getMyArgs(int argc, char* argv[], int& ll, int& lt, double& beta, double& t, double& U, double& mu,
               int& nwrap, int& nwarm, int& nbin, int& nsweep, int& nBetweenBins, std::string& filename_eqtime,
               std::string& filename_dynamic, bool& bool_Append, bool& bool_eqtime, bool& bool_dynamic) {

    std::stringstream sstream;
    std::vector<std::string> myArgs_temp(argc);

    // read parameters from command line
    for (int arg = 1; arg < argc; ++arg) {
        sstream << argv[arg];
        sstream >> myArgs_temp[arg];
        sstream.clear();
    }

    // delete empty strings
    std::vector<std::string> myArgs;
    for (const auto &arg : myArgs_temp) {
        if (!arg.empty()) {
            myArgs.push_back(arg);
        }
    }

    // make optionsMap
    std::map<std::string, std::string> optionsMap;
    while (!myArgs.empty()) {
        if (myArgs[0][0] == '-') {
            std::string str_temp = myArgs[0].substr(1, myArgs[0].size());
            myArgs[0] = str_temp;

            std::string key = myArgs[0];
            std::string value;

            // std::vector<std::string> table(2);
            // table[0] = myArgs[0];

            std::reverse(myArgs.begin(), myArgs.end());
            myArgs.pop_back();
            std::reverse(myArgs.begin(), myArgs.end());

            if (myArgs.empty()) {
                // table[1] = myArgs[0];
                value = "";
                // myTable.push_back(table);
                optionsMap[key] = value;
                break;
            } else {
                if (myArgs[0][0] != '-') {
                    value = myArgs[0];
                    std::reverse(myArgs.begin(), myArgs.end());
                    myArgs.pop_back();
                    std::reverse(myArgs.begin(), myArgs.end());
                } else { value = ""; }
                optionsMap[key] = value;
            }
        } else {
            std::reverse(myArgs.begin(), myArgs.end());
            myArgs.pop_back();
            std::reverse(myArgs.begin(), myArgs.end());
        }
    }

    // read parameters from optionsMap
    // help message first
    for (auto & iter : optionsMap) {
        if (option2int(iter.first) == 16 || option2int(iter.first) == 17)
            usage(argv[0]);
    }

    for (auto & iter : optionsMap) {
        std::string key = iter.first;
        std::string value = iter.second;

        if (!value.empty()) {
            switch (option2int(key)) {
                case 0: {
                    ll = str2int(value);
                    nwarm = ceil(4 * ll * ll * beta);
                    break;
                }
                case 1: {
                    lt = str2int(value);
                    break;
                }
                case 2: {
                    beta = str2double(value);
                    nwarm = ceil(4 * ll * ll * beta);
                    break;
                }
                case 3: {
                    t = str2double(value);
                    break;
                }
                case 4: {
                    U = str2double(value);
                    break;
                }
                case 5: {
                    mu = str2double(value);
                    break;
                }
                case 6: {
                    nwrap = str2int(value);
                    break;
                }
                case 8: {
                    nbin = str2int(value);
                    break;
                }
                case 9: {
                    nsweep = str2int(value);
                    break;
                }
                case 10: {
                    nBetweenBins = str2int(value);
                    break;
                }
                case 11: {
                    if (value == "true") { bool_Append = true; }
                    if (value == "false") { bool_Append = false; }
                    break;
                }
                case 12: {
                    if (value == "true") { bool_eqtime = true; }
                    if (value == "false") { bool_eqtime = false; }
                    break;
                }
                case 13: {
                    if (value == "true") { bool_dynamic = true; }
                    if (value == "false") { bool_dynamic = false; }
                    break;
                }
                case 14: {
                    filename_eqtime = value;
                    break;
                }
                case 15: {
                    filename_dynamic = value;
                    break;
                }
            }
        }
    }

    // read nwarm information at last
    for (auto & iter : optionsMap) {
        if (!iter.second.empty()) {
            if (option2int(iter.first) == 7)
                nwarm = str2int(iter.second);
        }
    }
}

#endif //DQMC_HUBBARD_OPTIONS_H
