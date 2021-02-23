#ifndef DQMC_HUBBARD_DETQMC_H
#define DQMC_HUBBARD_DETQMC_H
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

    double MomentumDist = 0.0;
    double localSpinCorr = 0.0;

    vecXd q = vecXd::Zero(2);

    time_t  begin_t, end_t;  // time cost of one single measuring process


public:

    detQMC() = default;

    void set_Model_Params(int ll, int lt, double beta, double t, double Uint, double mu, int nwrap=8);

    void set_MC_Params(int nwarm, int nsweep);

    void set_Momentum_q(double qx, double qy);

    void runQMC(bool bool_warm, bool bool_display_process);

    void printStats();

    void outputStats(const std::string& filename, bool bool_Append);


private:

    void sweep_BackAndForth(bool bool_measure);

    void measure_equal_time();

    /** double occupation: D = < n_up*n_dn > */
    double meas_DoubleOccu(const matXd& gu, const matXd& gd) const;

    /** single particle kinetic energy */
    double meas_KineticEnergy(const matXd& gu, const matXd& gd) const;

    /** momentum distribution of electrons: fourier transformation of real-space electron distribution */
    double meas_MomentumDist(const matXd& gu, const matXd& gd, const vecXd& p) const;

    /** local spin correlation: magnetization C(0,0) = < (n_up - n_dn)^2 > */
    double meas_localSpinCorr(const matXd& gu, const matXd& gd) const;

    /** magnetic struct factor: fourier transformation of real-space spin-spin correlation */
    double meas_StructFactor(const matXd& gu, const matXd& gd, const vecXd& p) const;


    /* todo */
    void measure_displaced_time();

    void clearStats();
};

#endif //DQMC_HUBBARD_DETQMC_H
