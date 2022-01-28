#ifndef DQMC_HUBBARD_MEASURE_GATHER_HPP
#define DQMC_HUBBARD_MEASURE_GATHER_HPP
#pragma once

/**
  *  This source file includes subroutines for efficiently 
  *  gathering measured observable data over different processors.
  *  Directly send/receive vectors of serialized data types or class objects.
  */ 


#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include "eigen_boost_serialization.hpp"
#include "observable.h"

namespace Measure {
    
    class GatherMPI {
        public:

        template<typename DataType>
        Observable<DataType> gather_observable(const boost::mpi::communicator &world, const Observable<DataType> &obs) {
            const int master = 0;
            const int rank = world.rank();

            // collect data from all processors
            if (rank == master) {
                std::vector<DataType> collected_data;
                std::vector<std::vector<DataType>> tmp_data(world.size()-1);
                std::vector<boost::mpi::request> recvs;

                // master processor
                collected_data.insert(collected_data.end(), obs.bin_data().begin(), obs.bin_data().end());

                // receive messages from other processors
                // pass vectors of serialized objects directly
                for (int proc = 1; proc < world.size(); ++proc) {
                    recvs.push_back(world.irecv(proc, proc, tmp_data[proc-1]));
                }
                boost::mpi::wait_all(recvs.begin(), recvs.end());
                for (auto data : tmp_data) {
                    collected_data.insert(collected_data.end(), data.begin(), data.end());
                }

                Observable<DataType> collected_obs;
                collected_obs.set_size_of_bin(collected_data.size());
                collected_obs.set_zero_element(obs.zero_element());
                collected_obs.set_observable_name(obs.name());
                // no need to add methods
                collected_obs.allocate();
                collected_obs.bin_data() = collected_data;

                // analysis
                collected_obs.analyse();
                return collected_obs;
            }
            else {
                std::vector<boost::mpi::request> sends;
                sends.push_back(world.isend(master, rank, obs.bin_data()));
                boost::mpi::wait_all(sends.begin(), sends.end());
                return obs;
            }
        }
        
        template<typename DataType>
        std::vector<Observable<DataType>> gather_observable(const boost::mpi::communicator &world, 
        const std::vector<Observable<DataType>> &obs_vec) {
            std::vector<Observable<DataType>> collected_obs_vec;
            for (auto obs : obs_vec) {
                collected_obs_vec.push_back(gather(world, obs));
            }
            return collected_obs_vec;
        }

        void gather_from_all_processors(const boost::mpi::communicator &world, Measure &measure) {
            const int master = 0;
            const int rank = world.rank();

            if (measure.is_eqtime_measure()) {
                *measure._container._sign_eqtime = gather_observable(world, *measure._container._sign_eqtime);
            }
            if (measure.is_dynamic_measure()) {
                *measure._container._sign_dynamic = gather_observable(world, *measure._container._sign_dynamic);
            }

            if (measure._container._obs_double) {
                for (auto &obs : *measure._container._obs_double) {
                    obs = gather_observable(world, obs);
                }
            }
            if (measure._container._obs_vector) {
                for (auto &obs : *measure._container._obs_vector) {
                    obs = gather_observable(world, obs);
                }
            }
            if (measure._container._obs_matrix) {
                for (auto &obs : *measure._container._obs_matrix) {
                    obs = gather_observable(world, obs);
                }
            }

            if (rank == master) {
                measure.set_size_of_bin(measure.nbin()*world.size());
                measure.analyse_stats();
            }
        }

    };

} // namespace Measure

#endif //DQMC_HUBBARD_MEASURE_GATHER_HPP