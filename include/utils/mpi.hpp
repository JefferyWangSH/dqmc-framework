#ifndef UTILS_MPI_HPP
#define UTILS_MPI_HPP
#pragma once

/**
  *  This source file includes implementations of the special mpi::gather method,
  *  which is designed to collect Observable::Observable classes among a set of MPI processes.
  */

#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include "utils/eigen_boost_serialization.hpp"
#include "measure/observable.h"
#include "measure/measure_handler.h"


namespace Utils {

    
    // ------------------------------------------  Utils::MPI  -------------------------------------------------
    class MPI {

        public:

            // gather the bin data of an observable object from other processes.
            // the data are collected into the master process,
            // and the corresponding Observable class is changed in place
            template<typename ObsType>
            static void gather_observable( const boost::mpi::communicator &world, 
                                           Observable::Observable<ObsType>* obs )
            {   
                // rank of the current process
                const int master = 0;
                const int rank = world.rank();

                // collect data from all other processes
                if ( rank == master ) {
                    std::vector<std::vector<ObsType>> tmp_data(world.size()-1);
                    std::vector<boost::mpi::request> recvs;

                    // for master processor, receive messages from other processes
                    // sending vectors of serialized objects directly
                    for ( auto proc = 1; proc < world.size(); ++proc ) {
                        recvs.push_back( world.irecv(proc, proc, tmp_data[proc-1]) );
                    }
                    boost::mpi::wait_all( recvs.begin(), recvs.end() );
                    
                    // push back the received data
                    for ( const auto& data : tmp_data ) {
                        obs->bin_data().insert( obs->bin_data().end(), data.begin(), data.end() );
                    }

                    // resize the number of bins
                    obs->set_number_of_bins( obs->bin_data().size() );
                }
                else {
                    // for subject processes, send observable data to the master
                    std::vector<boost::mpi::request> sends;
                    sends.push_back( world.isend( master, rank, obs->bin_data() ) );
                    boost::mpi::wait_all( sends.begin(), sends.end() );
                }
            }


            // gather all observable objects in the measuring handler
            // note that the Utils::MPI class should be a friend class of Measure::MeasureHandler
            // to get access to the protected observable members
            static void mpi_gather( const boost::mpi::communicator &world, Measure::MeasureHandler& meas_handler )
            {   
                // scalar observables
                for ( auto& scalar_obs : meas_handler.m_eqtime_scalar_obs ) {
                    gather_observable( world, scalar_obs.get() );
                }
                for ( auto& scalar_obs : meas_handler.m_dynamic_scalar_obs ) {
                    gather_observable( world, scalar_obs.get() );
                }

                // vector observables
                for ( auto& vector_obs : meas_handler.m_eqtime_vector_obs ) {
                    gather_observable( world, vector_obs.get() );
                }
                for ( auto& vector_obs : meas_handler.m_dynamic_vector_obs ) {
                    gather_observable( world, vector_obs.get() );
                }

                // matrix observables
                for ( auto& matrix_obs : meas_handler.m_eqtime_matrix_obs ) {
                    gather_observable( world, matrix_obs.get() );
                }
                for ( auto& matrix_obs : meas_handler.m_dynamic_matrix_obs ) {
                    gather_observable( world, matrix_obs.get() );
                }

                // statistics of the configuration sign
                if ( meas_handler.m_equaltime_sign ) {
                    gather_observable( world, meas_handler.m_equaltime_sign.get() );
                }
                if ( meas_handler.m_dynamic_sign ) {
                    gather_observable( world, meas_handler.m_dynamic_sign.get() );
                }

                // reset the number of bins
                meas_handler.m_bin_num *= world.size();
            }


    };


} // namespace Utils


#endif // UTILS_MPI_HPP
