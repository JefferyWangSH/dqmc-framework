#ifndef DQMC_HUBBARD_MEASURE_GATHER_HPP
#define DQMC_HUBBARD_MEASURE_GATHER_HPP
#pragma once

/**
  *  This source file includes the subroutine to efficiently gather measured observable data over different processors
  *  without realization of serialization of user-defined data types.
  *  TODO: transmit STL or user-defined data type directly ( need further tests ).
  */ 


#include <boost/mpi.hpp>
#include "measure_data.h"

namespace Measure {
    Measure::MeasureData gather(const boost::mpi::communicator &world, const Measure::MeasureData &obs) {
        const int master = 0;
        const int rank = world.rank();

        // collect data from all processors
        if (rank == master) {
            std::vector<double> bin_data(world.size() * obs.size_of_bin());
            std::vector<boost::mpi::request> recvs;

            // loop for bins and processors
            for (int proc = 0; proc < world.size(); ++proc) {
                for (int bin = 0; bin < obs.size_of_bin(); ++bin) {
                    if (proc == master) {
                        bin_data[bin] = obs.bin_data()[bin];
                    }
                    else {
                        recvs.push_back(world.irecv(proc, bin, bin_data[bin + proc * obs.size_of_bin()]));
                    }
                }
            }
            boost::mpi::wait_all(recvs.begin(), recvs.end());

            Measure::MeasureData gathered_obs(bin_data.size());
            gathered_obs.bin_data() = Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(bin_data.data(), 1, bin_data.size());

            // analysis
            gathered_obs.analyse();
            return gathered_obs;
        }
        else {
            std::vector<boost::mpi::request> sends;
            for (int bin = 0; bin < obs.size_of_bin(); ++bin) {
                sends.push_back(world.isend(master, bin, obs.bin_data()[bin]));
            }
            boost::mpi::wait_all(sends.begin(), sends.end());
            return obs;
        }
    }
    
} // namespace Measure

#endif //DQMC_HUBBARD_MEASURE_GATHER_HPP