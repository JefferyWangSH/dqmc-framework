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
#include "measure_data.h"

namespace Measure {

    template<typename DataStructure>
    Measure::MeasureData<DataStructure> gather(
        const boost::mpi::communicator &world, const Measure::MeasureData<DataStructure> &obs) {
        const int master = 0;
        const int rank = world.rank();

        // collect data from all processors
        if (rank == master) {
            std::vector<DataStructure> collected_data;
            std::vector<std::vector<DataStructure>> tmp_data(world.size()-1);
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

            Measure::MeasureData<DataStructure> collected_obs;
            collected_obs.set_size_of_bin(collected_data.size());
            collected_obs.set_zero_element(obs.zero_element());
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
    

    template<typename DataStructure>
    std::vector<Measure::MeasureData<DataStructure>> gather(
        const boost::mpi::communicator &world, 
        const std::vector<Measure::MeasureData<DataStructure>> &obs_vec) {
        
        std::vector<Measure::MeasureData<DataStructure>> collected_obs_vec;
        for (auto obs : obs_vec) {
            collected_obs_vec.push_back(Measure::gather(world, obs));
        }
        return collected_obs_vec;
    }

} // namespace Measure

#endif //DQMC_HUBBARD_MEASURE_GATHER_HPP