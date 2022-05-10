#ifndef OBSERABLE_HANDLER_H
#define OBSERABLE_HANDLER_H
#pragma once

/**
  *  This header file defines the Observable::ObservableHandler class 
  *  to handle with the list of observables during dqmc simulations.
  *  Automatic classification of input observables and the allocation of memory are included.
  */

#include <map>
#include <memory>
#include "measure/observable.h"


namespace Observable {

    using ObsName = std::string;
    using ObsNameList = std::vector<std::string>;
    using EqtimeObsNameList = std::vector<std::string>;
    using DynamicObsNameList = std::vector<std::string>;
    

    // ----------------------------- Handler class Observable::ObservableHandler ------------------------------
    class ObservableHandler {
        protected:

            using ObsMap = std::map<std::string, std::shared_ptr<ObservableBase>>;

            using ptrBaseObs = std::shared_ptr<ObservableBase>;
            using ptrScalarObs = std::shared_ptr<ScalarObs>;
            using ptrVectorObs = std::shared_ptr<VectorObs>;
            using ptrMatrixObs = std::shared_ptr<MatrixObs>;
            
            using EqtimeScalarObs = std::vector<std::shared_ptr<ScalarObs>>;
            using EqtimeVectorObs = std::vector<std::shared_ptr<VectorObs>>;
            using EqtimeMatrixObs = std::vector<std::shared_ptr<MatrixObs>>;
            using DynamicScalarObs = std::vector<std::shared_ptr<ScalarObs>>;
            using DynamicVectorObs = std::vector<std::shared_ptr<VectorObs>>;
            using DynamicMatrixObs = std::vector<std::shared_ptr<MatrixObs>>;

            // map of observables for quick reference
            // only for finding or searching cetain observable, and frequently calls should be avoided
            ObsMap m_obs_map{};

            EqtimeScalarObs m_eqtime_scalar_obs{};
            EqtimeVectorObs m_eqtime_vector_obs{};
            EqtimeMatrixObs m_eqtime_matrix_obs{};

            DynamicScalarObs m_dynamic_scalar_obs{};
            DynamicVectorObs m_dynamic_vector_obs{};
            DynamicMatrixObs m_dynamic_matrix_obs{};

            ptrScalarObs m_equaltime_sign{};
            ptrScalarObs m_dynamic_sign{};

            // list of supported physical observables
            EqtimeObsNameList m_eqtime_obs_name =   {
                                                        "filling_number", 
                                                        "double_occupancy", 
                                                        "kinetic_energy", 
                                                        "momentum_distribution", 
                                                        "local_spin_corr", 
                                                        "spin_density_structure_factor", 
                                                        "charge_density_structure_factor", 
                                                        "s_wave_pairing_corr",
                                                    };

            DynamicObsNameList m_dynamic_obs_name = {
                                                        "greens_functions", 
                                                        "density_of_states", 
                                                        "superfluid_stiffness", 
                                                    };

        public:

            ObservableHandler() = default;

            // check if certain observable exists
            bool find(const ObsName& obs_name);

            // return certain type of the observable class
            template<typename ObsType> const ObsType find(const ObsName& obs_name);

            // initialize the handler
            void initial(const ObsNameList& obs_list);


        private:
            
            // check if certain observable is of eqtime/dynamic type
            bool is_eqtime(const ObsName& obs_name) const;
            bool is_dynamic(const ObsName& obs_name) const;

            // check the validity of the input list of observables
            bool check_validity(const ObsNameList& obs_list) const;

            // deallocate memory
            void deallocate();

    };


    // implementation of template member function
    template<typename ObsType> 
    const ObsType ObservableHandler::find(const ObsName& obs_name) {
        if ( this->find(obs_name) ) {
            auto ptrObs = std::dynamic_pointer_cast<ObsType>(this->m_obs_map[obs_name]);
            if ( ptrObs ) { return *ptrObs; }
        }
        return ObsType();
    }


} // namespace Observable


#endif // OBSERVABLE_HANDLER_H