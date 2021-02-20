#ifndef HUBBARD_V1_3_DETQMC_H
#define HUBBARD_V1_3_DETQMC_H

#pragma once
#include "hubbard.h"

class detQMC
{
private:
    Hubbard hubb;
    int nwrap{8}, nwarm{300}, nsweep{200};

    // for measuring
    int nn = 0;
    double DoubleOccu = 0.0;
    double KineticEnergy = 0.0;
    double StructFactor = 0.0;

    time_t  begin_t, end_t;  // time cost of one single measuring process


public:

    detQMC() = default;

    void set_Model_Params(int ll, int lt, double beta, double t, double Uint, double mu, int nwrap=8);

    void set_MC_Params(int nwarm, int nsweep);

    void runQMC(bool bool_display_process);

    void printStats();

    void outputStats(const std::string& filename, bool bool_Append);


private:

    void sweep_BackAndForth(bool bool_measure);

    void measure();

    void clearStats();
};

#endif //HUBBARD_V1_3_DETQMC_H
