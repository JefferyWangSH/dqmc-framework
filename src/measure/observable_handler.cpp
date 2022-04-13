#include "measure/observable_handler.h"
#include "measure/measure_methods.h"
#include <iostream>


namespace Observable {


    bool ObservableHandler::is_eqtime(const ObsName& obs_name) const {
        return ( std::find( this->m_eqtime_obs_name.begin(), 
                            this->m_eqtime_obs_name.end(), 
                            obs_name  ) 
                != this->m_eqtime_obs_name.end() );
    }

    bool ObservableHandler::is_dynamic(const ObsName& obs_name) const {
        return ( std::find( this->m_dynamic_obs_name.begin(), 
                            this->m_dynamic_obs_name.end(), 
                            obs_name  ) 
                != this->m_dynamic_obs_name.end() );
    }

    bool ObservableHandler::find(const ObsName& obs_name) {
        return ( this->m_obs_map.find(obs_name) != this->m_obs_map.end() );
    }

    const ScalarObs ObservableHandler::find_scalar(const ObsName& obs_name) {
        if ( this->find(obs_name) ) {
            auto ptrObs = std::dynamic_pointer_cast<ScalarObs>(this->m_obs_map[obs_name]);
            if ( ptrObs ) { return *ptrObs; }
        }
        return ScalarObs();
    }
    
    const VectorObs ObservableHandler::find_vector(const ObsName& obs_name) {
        if ( this->find(obs_name) ) {
            auto ptrObs = std::dynamic_pointer_cast<VectorObs>(this->m_obs_map[obs_name]);
            if ( ptrObs ) { return *ptrObs; }
        }
        return VectorObs();
    }

    const MatrixObs ObservableHandler::find_matrix(const ObsName& obs_name) {
        if ( this->find(obs_name) ) {
            auto ptrObs = std::dynamic_pointer_cast<MatrixObs>(this->m_obs_map[obs_name]);
            if ( ptrObs ) { return *ptrObs; }
        }
        return MatrixObs();
    }


    bool ObservableHandler::check_validity(const ObsNameList& obs_list) const {
        // preprocessing
        // remove redundant input
        ObsNameList tmp_list = obs_list;
        std::sort(tmp_list.begin(), tmp_list.end());
        tmp_list.erase(unique(tmp_list.begin(),tmp_list.end()), tmp_list.end());

        // check the validity of the input
        for (const auto& obs : tmp_list) {
            if ( !this->is_eqtime(obs) && !this->is_dynamic(obs) ) {
                return false;
            }
        }
        // otherwise return true
        return true;
    }


    void ObservableHandler::initial(const ObsNameList& obs_list) {
        
        // check the validity of the input
        if ( !this->check_validity(obs_list) ) {
            // unsupported observables found, throw errors
            std::cerr << " Unsupported observable type from the input. " << std::endl;
            exit(1);
        }

        for (const auto& obs_name : obs_list) {
            // allocate memory for observables
            // caution that to this stage only properties like name or method is assigned.
            // other info, such as m_size_of_bin and dimensional of m_zero_elem,
            // is kept unassigned until MeasureHandler or Lattice class is specialized.


            // ------------------------- Equal-time Observables --------------------------

            // ---------------------------- Filling number -------------------------------
            if ( obs_name == "filling_number" ) {
                ptrScalarObs filling_number = std::make_shared<ScalarObs>();
                filling_number->set_observable_name(obs_name);
                filling_number->add_method(Measure::Methods::measure_filling_number);
                this->m_eqtime_scalar_obs.emplace_back(filling_number);

                // fill in the map
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>(filling_number);
            }

            // --------------------------- Double occupancy ------------------------------
            if ( obs_name == "double_occupancy" ) {
                ptrScalarObs double_occupancy = std::make_shared<ScalarObs>();
                double_occupancy->set_observable_name(obs_name);
                double_occupancy->add_method(Measure::Methods::measure_double_occupancy);
                this->m_eqtime_scalar_obs.emplace_back(double_occupancy);
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>(double_occupancy);
            }

            // --------------------------- Kinetic energy --------------------------------
            if ( obs_name == "kinetic_energy" ) {
                ptrScalarObs kinetic_energy = std::make_shared<ScalarObs>();
                kinetic_energy->set_observable_name(obs_name);
                kinetic_energy->add_method(Measure::Methods::measure_kinetic_energy);
                this->m_eqtime_scalar_obs.emplace_back(kinetic_energy);
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>(kinetic_energy);
            }

            // ------------------------ Momentum distribution ----------------------------
            if ( obs_name == "momentum_distribution" ) {
                ptrScalarObs momentum_distribution = std::make_shared<ScalarObs>();
                momentum_distribution->set_observable_name(obs_name);
                momentum_distribution->add_method(Measure::Methods::measure_momentum_distribution);
                this->m_eqtime_scalar_obs.emplace_back(momentum_distribution);
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>(momentum_distribution);
            }

            // ------------------------ Local spin correlations --------------------------
            if ( obs_name == "local_spin_corr" ) {
                ptrScalarObs local_spin_corr = std::make_shared<ScalarObs>();
                local_spin_corr->set_observable_name(obs_name);
                local_spin_corr->add_method(Measure::Methods::measure_local_spin_corr);
                this->m_eqtime_scalar_obs.emplace_back(local_spin_corr);
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>(local_spin_corr);
            }

            // --------------- Spin density wave (SDW) structure factor ------------------
            if ( obs_name == "spin_density_structure_factor" ) {
                ptrScalarObs sdw_factor = std::make_shared<ScalarObs>();
                sdw_factor->set_observable_name(obs_name);
                sdw_factor->add_method(Measure::Methods::measure_spin_density_structure_factor);
                this->m_eqtime_scalar_obs.emplace_back(sdw_factor);
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>(sdw_factor);
            }

            // ------------- Charge density wave (CDW) structure factor ------------------
            if ( obs_name == "charge_density_structure_factor" ) {
                ptrScalarObs cdw_factor = std::make_shared<ScalarObs>();
                cdw_factor->set_observable_name(obs_name);
                cdw_factor->add_method(Measure::Methods::measure_charge_density_structure_factor);
                this->m_eqtime_scalar_obs.emplace_back(cdw_factor);
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>(cdw_factor);
            }

            // ----------- S wave pairing correlations of superconductivity --------------
            if ( obs_name == "s_wave_pairing_corr" ) {
                ptrScalarObs s_wave_pairing_corr = std::make_shared<ScalarObs>();
                s_wave_pairing_corr->set_observable_name(obs_name);
                s_wave_pairing_corr->add_method(Measure::Methods::measure_s_wave_pairing_corr);
                this->m_eqtime_scalar_obs.emplace_back(s_wave_pairing_corr);
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>(s_wave_pairing_corr);
            }

            // adding new methods here


            // --------------------- Time-displaced Observables --------------------------

            // --------------------- Greens functions G(k, tau) --------------------------
            if ( obs_name == "greens_functions" ) {
                ptrMatrixObs greens_functions = std::make_shared<MatrixObs>();
                greens_functions->set_observable_name(obs_name);
                greens_functions->add_method(Measure::Methods::measure_greens_functions);
                this->m_dynamic_matrix_obs.emplace_back(greens_functions);
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>(greens_functions);
            }

            // --------------------- Density of states D(tau) ----------------------------
            if ( obs_name == "density_of_states" ) {
                ptrVectorObs density_of_states = std::make_shared<VectorObs>();
                density_of_states->set_observable_name(obs_name);
                density_of_states->add_method(Measure::Methods::measure_density_of_states);
                this->m_dynamic_vector_obs.emplace_back(density_of_states);
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>(density_of_states);
            }


            // ---------------------- Superfluid stiffness -------------------------------
            if ( obs_name == "superfluid_stiffness" ) {
                ptrScalarObs superfluid_stiffness = std::make_shared<ScalarObs>();
                superfluid_stiffness->set_observable_name(obs_name);
                superfluid_stiffness->add_method(Measure::Methods::measure_superfluid_stiffness);
                this->m_dynamic_scalar_obs.emplace_back(superfluid_stiffness);
                this->m_obs_map[obs_name] = std::static_pointer_cast<ObservableBase>(superfluid_stiffness);
            }

            // add new methods here
        
        }

        // adding measurements of configuration signs manually to keep track of the sign problem
        if (    !this->m_eqtime_scalar_obs.empty() 
             || !this->m_eqtime_vector_obs.empty() 
             || !this->m_eqtime_matrix_obs.empty() ) {
            ptrScalarObs eqtime_sign = std::make_shared<ScalarObs>();
            eqtime_sign->set_observable_name("eqtime_sign");
            eqtime_sign->add_method(Measure::Methods::measure_eqtime_config_sign);
            this->m_eqtime_scalar_obs.emplace_back(eqtime_sign);
            this->m_obs_map["eqtime_sign"] = std::static_pointer_cast<ObservableBase>(eqtime_sign);
        }

        if (    !this->m_dynamic_scalar_obs.empty() 
             || !this->m_dynamic_vector_obs.empty() 
             || !this->m_dynamic_matrix_obs.empty() ) {
            ptrScalarObs dynamic_sign = std::make_shared<ScalarObs>();
            dynamic_sign->set_observable_name("dynamic_sign");
            dynamic_sign->add_method(Measure::Methods::measure_dynamic_config_sign);
            this->m_dynamic_scalar_obs.emplace_back(dynamic_sign);
            this->m_obs_map["dynamic_sign"] = std::static_pointer_cast<ObservableBase>(dynamic_sign);
        }
    }


} // namespace Observable