#ifndef DQMC_HUBBARD_MEASURE_CONTAINER_H
#define DQMC_HUBBARD_MEASURE_CONTAINER_H
#pragma once

/**
  *  This header file includes Measure::Container class 
  *  for collection of observables to be measured.
  *  Automatic classification of input observables and allocating process supported.
  */

#include <memory>
#include "observable.h"


namespace Measure {

    class Container {
    private:
        using StringList = std::vector<std::string>;
        using DoubleObservableList = std::vector<Observable<double>>;
        using VectorObservableList = std::vector<Observable<Eigen::VectorXd>>;
        using MatrixObservableList = std::vector<Observable<Eigen::MatrixXd>>;

        std::unique_ptr<StringList> _obs_list{};
        std::unique_ptr<StringList> _obs_list_eqtime{};
        std::unique_ptr<StringList> _obs_list_dynamic{};
        std::unique_ptr<DoubleObservableList> _obs_double{};
        std::unique_ptr<VectorObservableList> _obs_vector{};
        std::unique_ptr<MatrixObservableList> _obs_matrix{};
        std::unique_ptr<Observable<double>> _sign_eqtime{}, _sign_dynamic{};

        // define all supported physical observables for measurements
        StringList _supported_obs_eqtime = {"filling_number", 
                                            "double_occupancy", 
                                            "kinetic_energy", 
                                            "momentum_distribution", 
                                            "local_spin_corr", 
                                            "spin_density_structure_factor", 
                                            "charge_density_structure_factor", 
                                            "s_wave_pairing_corr",
                                            };

        StringList _supported_obs_dynamic = {"greens_functions", 
                                            "density_of_states", 
                                            "superfluid_stiffness", 
                                            };

        friend class Measure;
        friend class GatherMPI;
    
    public:
        // (de)constructions 
        Container() = default;
        ~Container() = default;

        // interface for the output of observables
        StringList list() const;
        StringList eqtime_list() const;
        StringList dynamic_list() const;
        DoubleObservableList obs_list_double() const;
        VectorObservableList obs_list_vector() const;
        MatrixObservableList obs_list_matrix() const;

        bool is_eqtime_measure() const;
        bool is_dynamic_measure() const;
        bool is_eqtime_obs(std::string obs) const;
        bool is_dynamic_obs(std::string obs) const;

        // read list of observables which are to be measured
        void read_list(const StringList &obs_list);

        // segment input observables into two classes, equal-time and time-displaced (dynamic)
        void filter();

        // allocate Observable class according to the input list
        void allocate();

    };

} // namespace Measure

#endif //DQMC_HUBBARD_MEASURE_CONTAINER_H