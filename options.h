#ifndef HUBBARD_V1_3_OPTIONS_H
#define HUBBARD_V1_3_OPTIONS_H

#include <map>

int option2int(const std::string& option) {
    if ( option == "ll") { return 0;}
    if ( option == "lt") { return 1;}
    if ( option == "beta") { return 2;}
    if ( option == "t") { return 3;}
    if ( option == "u") { return 4;}
    if ( option == "mu") { return 5;}
    if ( option == "nwrap") { return 6;}
    if ( option == "nwarm") { return 7;}
    if ( option == "nsweep") { return 8;}
    if ( option == "o") { return 9;}
    if ( option == "app") { return 10;}
    if ( option == "h") { return 11;}
    if ( option == "?") { return 12;}
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
              << "   [ -ll  4 ]         (spatial) size of lattice, default: 4\n"
              << "   [ -lt  80 ]        (imaginary time) size of lattice, default: 80\n"
              << "   [ -beta  4.0 ]     inverse temperature, default: 4.0\n"
              << "   [ -t  1.0 ]        hopping strength: 1.0\n"
              << "   [ -u  4.0 ]        interaction strength, default: 4.0\n"
              << "   [ -mu  0.0 ]       chemical potential, default: 0.0\n"
              << "   [ -nwrap  10 ]     pace of stabilization process, default: 10\n"
              << "   [ -nwarm  256 ]    number of warmup sweeps, default: 256\n"
              << "   [ -nsweep  200 ]   number of measurement sweeps, default: 200\n"
              << "   [ -o  output.txt ] output filename, default: output.txt\n"
              << "   [ -app  true ]     outfile mode: app or trunc, default: true\n"
              << "   [ -h / -? ]        display this information\n"
              << std::endl;
    exit(1);
}

// function to read and reload parameters from the command line
void getMyArgs(int argc, char* argv[], int& ll, int& lt, double& beta, double& t, double& U, double& mu,
               int& nwrap, int& nwarm, int& nsweep, std::string& filename, bool& bool_Append) {

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
    for (auto iter = optionsMap.begin(); iter != optionsMap.end(); ++iter) {
        if (option2int(iter->first) == 11 || option2int(iter->first) == 12)
            usage(argv[0]);
    }

    for (auto iter = optionsMap.begin(); iter != optionsMap.end(); ++iter) {
        std::string key = iter->first;
        std::string value = iter->second;

        if (!value.empty()) {
            switch (option2int(key)) {
                case 0: {
                    ll = str2int(value);
                    nwarm = 4 * ll * ll * beta;
                    break;
                }
                case 1: {
                    lt = str2int(value);
                    break;
                }
                case 2: {
                    beta = str2double(value);
                    nwarm = 4 * ll * ll * beta;
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
                    nsweep = str2int(value);
                    break;
                }
                case 9: {
                    filename = value;
                    break;
                }
                case 10: {
                    if (value == "true") { bool_Append = true; }
                    if (value == "false") { bool_Append = false; }
                    break;
                }
            }
        }
    }

    // read nwarm statistics at last
    for (auto iter = optionsMap.begin(); iter != optionsMap.end(); ++iter) {
        if (!iter->second.empty()) {
            if (option2int(iter->first) == 7)
                nwarm = str2int(iter->second);
        }
    }
}

#endif //HUBBARD_V1_3_OPTIONS_H
