#ifndef DQMC_IO_H
#define DQMC_IO_H
#pragma once

/**
  *  This header file defines QuantumMonteCarlo::DqmcIO class 
  *  containing the basic IO interfaces for the input/output of the dqmc data
  */


#include <fstream>
#include <string>
#include <boost/format.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "dqmc.h"
#include "dqmc_walker.h"

#include "model/model_base.h"
#include "model/repulsive_hubbard.h"
#include "model/attractive_hubbard.h"

#include "lattice/lattice_base.h"
#include "lattice/square.h"
#include "lattice/cubic.h"
#include "lattice/honeycomb.h"

#include "checkerboard/checkerboard_base.h"
#include "measure/measure_handler.h"
#include "measure/observable.h"


namespace QuantumMonteCarlo {

    using ModelBase = Model::ModelBase;
    using LatticeBase = Lattice::LatticeBase;
    using MeasureHandler = Measure::MeasureHandler;
    using CheckerBoardBasePtr = std::unique_ptr<CheckerBoard::CheckerBoardBase>;


    // ---------------------------- IO Interface class QuantumMonteCarlo::DqmcIO ------------------------------
    class DqmcIO {

        public:

            // output the information of dqmc initialization,
            // including initialization status and simulation parameters.
            // the behavior of this function depends on specific model and lattice types.
            template<typename StreamType>
            static void output_init_info              ( StreamType& ostream, 
                                                        const ModelBase& model, 
                                                        const LatticeBase& lattice,
                                                        const DqmcWalker& walker, 
                                                        const MeasureHandler& meas_handler,
                                                        const CheckerBoardBasePtr& checkerboard );
            
            // output the ending information of the simulation,
            // including time cost and wrapping errors
            template<typename StreamType>
            static void output_ending_info            ( StreamType& ostream, const DqmcWalker& walker );

            // output the mean value and error bar of one specific observable
            template<typename StreamType, typename ObsType>
            static void output_observable             ( StreamType& ostream, const Observable::Observable<ObsType>& obs );

            // output the bin data of one specific observable
            template<typename StreamType, typename ObsType>
            static void output_observable_in_bins     ( StreamType& ostream, const Observable::Observable<ObsType>& obs );
            
            // output list of inequivalent momentum points ( k stars )
            template<typename StreamType>
            static void output_k_stars                ( StreamType& ostream, const LatticeBase& lattice );

            // output imgainary-time grids
            template<typename StreamType>
            static void output_imaginary_time_grids   ( StreamType& ostream, const DqmcWalker& walker );

            // output the current configuration the bosonic fields,
            // depending on specific model type
            template<typename StreamType>
            static void output_bosonic_fields         ( StreamType& ostream, const ModelBase& model );
            
            // read the configuration of the bosonic fields from input file
            // depending on specific model type
            static void read_bosonic_fields_from_file ( const std::string& filename, ModelBase& model );

    };





    // -------------------------------------------------------------------------------------------------------
    //                                Implementation of template functions
    // -------------------------------------------------------------------------------------------------------


    template<typename StreamType>
    void DqmcIO::output_init_info( StreamType& ostream, 
                                   const ModelBase& model, 
                                   const LatticeBase& lattice, 
                                   const DqmcWalker& walker,
                                   const MeasureHandler& meas_handler,
                                   const CheckerBoardBasePtr& checkerboard )
    {
        if ( !ostream ) {
            std::cerr << "QuantumMonteCarlo::DqmcIO::output_init_info(): "
                      << "the ostream failed to work, please check the input." << std::endl;
            exit(1);
        }
        else {
            
            // output formats
            boost::format fmt_param_str("%| 28s|%| 7s|%| 24s|\n");
            boost::format fmt_param_int("%| 28s|%| 7s|%| 24d|\n");
            boost::format fmt_param_double("%| 28s|%| 7s|%| 24.3f|\n");
            const std::string_view joiner = "->";
            auto bool2str = [](bool b) {if (b) return "True"; else return "False";};
            
            // -------------------------------------------------------------------------------------------
            //                             Output current date and time
            // -------------------------------------------------------------------------------------------
            const auto current_time = boost::posix_time::second_clock::local_time();
            ostream << boost::format(" Current time: %s\n") % current_time << std::endl;


            // -------------------------------------------------------------------------------------------
            //                                 Output model information
            // -------------------------------------------------------------------------------------------

            // ------------------------------  Repulsive Hubbard model  ----------------------------------
            if ( const auto repulsive_hubbard = dynamic_cast<const Model::RepulsiveHubbard*>(&model);
                repulsive_hubbard != nullptr ) {
                ostream << " Model: Attractive Hubbard\n\n"
                        << fmt_param_double % "Hopping constant \'t\'" % joiner % repulsive_hubbard->HoppingT()
                        << fmt_param_double % "Onsite interaction \'U\'" % joiner % repulsive_hubbard->OnSiteU()
                        << fmt_param_double % "Checimcal potential \'mu\'" % joiner % repulsive_hubbard->ChemicalPotential()
                        << std::endl;
            }
            
            // ------------------------------  Attractive Hubbard model  ---------------------------------
            else if ( const auto attractive_hubbard = dynamic_cast<const Model::AttractiveHubbard*>(&model); 
                attractive_hubbard != nullptr ) {
                ostream << " Model: Attractive Hubbard\n\n"
                        << fmt_param_double % "Hopping constant \'t\'" % joiner % attractive_hubbard->HoppingT()
                        << fmt_param_double % "Onsite interaction \'U\'" % joiner % attractive_hubbard->OnSiteU()
                        << fmt_param_double % "Checimcal potential \'mu\'" % joiner % attractive_hubbard->ChemicalPotential()
                        << std::endl;
            }

            else {
                std::cerr << "QuantumMonteCarlo::DqmcIO::output_init_info(): "
                          << "undefined model type." << std::endl;
                exit(1);
            }


            // -------------------------------------------------------------------------------------------
            //                                Output lattice information
            // -------------------------------------------------------------------------------------------
            // output the momentum points by the way

            // ---------------------------------  2d Square lattice  -------------------------------------
            if ( const auto square_lattice = dynamic_cast<const Lattice::Square*>(&lattice);
                square_lattice != nullptr ) {
                boost::format fmt_cell("%d * %d");
                boost::format fmt_momentum("(%.2f, %.2f) pi");
                const int side_length = square_lattice->SideLength();
                const double px = (square_lattice->Index2Momentum(meas_handler.Momentum(), 0)/M_PI);
                const double py = (square_lattice->Index2Momentum(meas_handler.Momentum(), 1)/M_PI);

                ostream << " Lattice: Square lattice\n\n"
                        << fmt_param_str % "Size of cell" % joiner % ( fmt_cell % side_length % side_length )
                        << fmt_param_str % "Momentum point" % joiner % ( fmt_momentum % px % py )
                        << std::flush;
            }
            
            // // ----------------------------------  3d Cubic lattice  -------------------------------------
            // else if ( const auto cubic_lattice = dynamic_cast<const Lattice::Cubic*>(&lattice);
            //     cubic_lattice != nullptr ) {
            //     // todo
            // }

            // // --------------------------------  2d Honeycomb lattice  -----------------------------------
            // else if ( const auto honeycomb_lattice = dynamic_cast<const Lattice::Honeycomb*>(&lattice);
            //     honeycomb_lattice != nullptr ) {
            //     // todo
            // }

            else {
                std::cerr << "QuantumMonteCarlo::DqmcIO::output_init_info(): "
                          << "undefined lattice type." << std::endl;
                exit(1); 
            }


            // -------------------------------------------------------------------------------------------
            //                              Output CheckerBoard information
            // -------------------------------------------------------------------------------------------
            ostream << fmt_param_str % "Checkerboard breakups" % joiner % bool2str((bool)checkerboard)
                    << std::endl;


            // -------------------------------------------------------------------------------------------
            //                               Output MonteCarlo Params
            // -------------------------------------------------------------------------------------------
            ostream << " MonteCarlo Params:\n\n"
                    << fmt_param_double % "Inverse temperature" % joiner % walker.Beta()
                    << fmt_param_int % "Imaginary-time length" % joiner % walker.TimeSize()
                    << fmt_param_double % "Imaginary-time interval" % joiner % walker.TimeInterval()
                    << fmt_param_int % "Stabilization pace" % joiner % walker.StabilizationPace()
                    << std::endl;

            // -------------------------------------------------------------------------------------------
            //                                Output Measuring Params
            // -------------------------------------------------------------------------------------------
            ostream << " Measuring Params:\n\n"
                    << fmt_param_str % "Warm up" % joiner % bool2str(meas_handler.isWarmUp())
                    << fmt_param_str % "Equal-time measure" % joiner % bool2str(meas_handler.isEqualTime())
                    << fmt_param_str % "Dynamical measure" % joiner % bool2str(meas_handler.isDynamic())
                    << std::endl;
            
            ostream << fmt_param_int % "Sweeps for warmup" % joiner % meas_handler.WarmUpSweeps()
                    << fmt_param_int % "Number of bins" % joiner % meas_handler.BinsNum()
                    << fmt_param_int % "Sweeps per bin" % joiner % meas_handler.BinsSize()
                    << fmt_param_int % "Sweeps between bins" % joiner % meas_handler.SweepsBetweenBins()
                    << std::endl;
        }
    }


    template<typename StreamType>
    void DqmcIO::output_ending_info( StreamType& ostream, const DqmcWalker& walker )
    {
        if ( !ostream ) {
            std::cerr << "QuantumMonteCarlo::DqmcIO::output_ending_info(): "
                      << "the ostream failed to work, please check the input." << std::endl;
            exit(1);
        }
        else {
            // parse the duration
            const int day = std::floor((double)Dqmc::timer() / 86400000);
            const int hour = std::floor(((double)Dqmc::timer()/1000 - day * 86400) / 3600);
            const int minute = std::floor(((double)Dqmc::timer()/1000 - day * 86400 - hour * 3600) / 60);
            const double sec = (double)Dqmc::timer()/1000 - 86400 * day - 3600 * hour - 60 * minute;

            // output the time cost of simulation
            if ( day ) { ostream << boost::format("\n The simulation finished in %d d %d h %d m %.2f s. \n") % day % hour % minute % sec << std::endl; }
            else if ( hour ) { ostream << boost::format("\n The simulation finished in %d h %d m %.2f s. \n") % hour % minute % sec << std::endl; }
            else if ( minute ) { ostream << boost::format("\n The simulation finished in %d m %.2f s. \n") % minute % sec << std::endl; }
            else { ostream << boost::format("\n The simulation finished in %.2f s. \n") % sec << std::endl; }

            // output wrapping errors of the evaluations of Green's functions
            ostream << boost::format(" Maximum of the wrapping error: %.5e\n") % walker.WrapError() << std::endl;
        }
    }


    template<typename StreamType, typename ObsType>
    void DqmcIO::output_observable( StreamType& ostream, const Observable::Observable<ObsType>& obs )
    {
        if ( !ostream ) {
            std::cerr << "QuantumMonteCarlo::DqmcIO::output_observable(): "
                      << "the ostream failed to work, please check the input." << std::endl;
            exit(1);
        }
        else {
            // standard screen output
            if constexpr ( std::is_same_v<StreamType, std::ostream> ) {
                
                // for scalar observables
                if constexpr ( std::is_same_v<ObsType, Observable::ScalarType> ) {
                    boost::format fmt_scalar_obs("%| 28s|%| 7s|%| 20.12f|  pm  %.12f");
                    const std::string joiner = "->";
                    ostream << fmt_scalar_obs % obs.description() % joiner % obs.mean_value() % obs.error_bar() << std::endl;
                }

                // // todo: currently not used
                // // for vector observables
                // else if constexpr ( std::is_same_v<ObsType, Observable::VectorType> ) {
                    
                // }

                // // for matrix observables
                // else if constexpr ( std::is_same_v<ObsType, Observable::MatrixType> ) {
                    
                // }

                // other observable type, raising errors
                else {
                    std::cerr << "QuantumMonteCarlo::DqmcIO::output_observable(): "
                              << "undefined observable type." << std::endl;
                    exit(1);
                }
            }

            // standard file output
            else if constexpr ( std::is_same_v<StreamType, std::ofstream> ) {
                
                // for scalar observables
                if constexpr ( std::is_same_v<ObsType, Observable::ScalarType> ) {
                    // for specfic scalar observable, output the mean value, error bar and relative error in order. 
                    boost::format fmt_scalar_obs("%| 20.10f|%| 20.10f|%| 20.10f|");
                    ostream << fmt_scalar_obs % obs.mean_value() % obs.error_bar() % (obs.error_bar()/obs.mean_value()) << std::endl;
                }

                // for vector observables
                else if constexpr ( std::is_same_v<ObsType, Observable::VectorType> ) {
                    // output vector observable
                    boost::format fmt_size_info("%| 20d|");
                    boost::format fmt_vector_obs("%| 20d|%| 20.10f|%| 20.10f|%| 20.10f|");

                    const int size = obs.mean_value().size();
                    const auto relative_error = (obs.error_bar().array()/obs.mean_value().array()).matrix();
                    ostream << fmt_size_info % size << std::endl;
                    for ( auto i = 0; i < size; ++i ) {
                        // output the mean value, error bar and relative error in order. 
                        ostream << fmt_vector_obs % i % obs.mean_value()(i) % obs.error_bar()(i) % relative_error(i) << std::endl;
                    }
                }

                // for matrix observables
                else if constexpr ( std::is_same_v<ObsType, Observable::MatrixType> ) {
                    // output matrix observable 
                    boost::format fmt_size_info("%| 20d|%| 20d|");
                    boost::format fmt_matrix_obs("%| 20d|%| 20d|%| 20.10f|%| 20.10f|%| 20.10f|");

                    const int row = obs.mean_value().rows();
                    const int col = obs.mean_value().cols();
                    const auto relative_error = (obs.error_bar().array()/obs.mean_value().array()).matrix();
                    ostream << fmt_size_info % row % col << std::endl;
                    for ( auto i = 0; i < row; ++i ) {
                        for ( auto j = 0; j < col; ++j ) {
                            // output the mean value, error bar and relative error in order. 
                            ostream << fmt_matrix_obs % i % j % obs.mean_value()(i,j) % obs.error_bar()(i,j) % relative_error(i,j) << std::endl;
                        }
                    }
                }

                // other observable types, raising errors
                else {
                    std::cerr << "QuantumMonteCarlo::DqmcIO::output_observable(): "
                              << "undefined observable type." << std::endl;
                    exit(1);
                }
            }

            // others stream types, raising errors
            else {
                std::cerr << "QuantumMonteCarlo::DqmcIO::output_observable(): "
                          << "unsupported type of output stream." << std::endl;
                exit(1);
            }
        }
    }


    template<typename StreamType, typename ObsType>
    void DqmcIO::output_observable_in_bins( StreamType& ostream, const Observable::Observable<ObsType>& obs )
    {
        if ( !ostream ) {
            std::cerr << "QuantumMonteCarlo::DqmcIO::output_observable_in_bins(): "
                      << "the ostream failed to work, please check the input." << std::endl;
            exit(1);
        }
        else {
            // for scalar observables
            if constexpr ( std::is_same_v<ObsType, Observable::ScalarType> ) {
                // output bin data of scalar observable
                boost::format fmt_size_info("%| 20d|");
                boost::format fmt_scalar_obs("%| 20d|%| 20.10f|");

                const int number_of_bins = obs.bin_num();
                ostream << fmt_size_info % number_of_bins << std::endl;
                for ( auto bin = 0; bin < number_of_bins; ++bin ) {
                    ostream << fmt_scalar_obs % bin % obs.bin_data(bin) << std::endl;
                }
            }

            // for vector observables
            else if constexpr ( std::is_same_v<ObsType, Observable::VectorType> ) {
                // output bin data of vector observable
                boost::format fmt_size_info("%| 20d|%| 20d|");
                boost::format fmt_vector_obs("%| 20d|%| 20d|%| 20.10f|");

                const int number_of_bins = obs.bin_num();
                const int size = obs.mean_value().size();
                ostream << fmt_size_info % number_of_bins % size << std::endl;
                for ( auto bin = 0; bin < number_of_bins; ++bin ) {
                    for ( auto i = 0; i < size; ++i ) {
                        ostream << fmt_vector_obs % bin % i % obs.bin_data(bin)(i) << std::endl;
                    }
                }
            }

            // for matrix observables
            else if constexpr ( std::is_same_v<ObsType, Observable::MatrixType> ) {
                // output bin data of matrix observable
                boost::format fmt_size_info("%| 20d|%| 20d|%| 20d|");
                boost::format fmt_matrix_obs("%| 20d|%| 20d|%| 20d|%| 20.10f|");

                const int number_of_bins = obs.bin_num();
                const int row = obs.mean_value().rows();
                const int col = obs.mean_value().cols();
                ostream << fmt_size_info % number_of_bins % row % col << std::endl;
                for ( auto bin = 0; bin < number_of_bins; ++bin ) {
                    for ( auto i = 0; i < row; ++i ) {
                        for ( auto j = 0; j < col; ++j ) {
                            ostream << fmt_matrix_obs % bin % i % j % obs.bin_data(bin)(i,j) << std::endl;
                        }
                    }
                }
            }

            // other observable types, raising errors
            else {
                std::cerr << "QuantumMonteCarlo::DqmcIO::output_observable_in_bins(): "
                          << "undefined observable type." << std::endl;
                exit(1);
            }
        }
    }

    
    template<typename StreamType>
    void DqmcIO::output_k_stars( StreamType& ostream, const LatticeBase& lattice )
    {
        if ( !ostream ) {
            std::cerr << "QuantumMonteCarlo::DqmcIO::output_k_stars(): "
                      << "the ostream failed to work, please check the input." << std::endl;
            exit(1);
        }
        else {
            // output k stars list
            boost::format fmt_info("%| 20d|");
            boost::format fmt_kstars("%| 20.10f|");
            ostream << fmt_info % lattice.kStarsNum() << std::endl;
            // loop for inequivalent momentum points
            for ( auto i = 0; i < lattice.kStarsNum(); ++i ) {
                ostream << fmt_info % i;
                // loop for axises of the reciprocal lattice
                for ( auto axis = 0; axis < lattice.SpaceDim(); ++axis ) {
                    ostream << fmt_kstars % lattice.Index2Momentum(i, axis);
                }
                ostream << std::endl;
            }
        }
    }


    template<typename StreamType>
    void DqmcIO::output_imaginary_time_grids( StreamType& ostream, const DqmcWalker& walker )
    {
        if ( !ostream ) {
            std::cerr << "QuantumMonteCarlo::DqmcIO::output_imaginary_time_grids(): "
                      << "the ostream failed to work, please check the input." << std::endl;
            exit(1);
        }
        else {
            // output the imaginary-time grids
            boost::format fmt_tgrids_info("%| 20d|%| 20.5f|%| 20.5f|");
            boost::format fmt_tgrids("%| 20d|%| 20.10f|");
            ostream << fmt_tgrids_info % walker.TimeSize() % walker.Beta() % walker.TimeInterval() << std::endl;
            for ( auto t = 0; t < walker.TimeSize(); ++t ) {
                ostream << fmt_tgrids % t % ( t * walker.TimeInterval() ) << std::endl;
            }
        }
    }


    template<typename StreamType>
    void DqmcIO::output_bosonic_fields( StreamType& ostream, const ModelBase& model )
    {
        if ( !ostream ) {
            std::cerr << "QuantumMonteCarlo::DqmcIO::output_bosonic_fields(): "
                      << "the ostream failed to work, please check the input." << std::endl;
            exit(1);
        }
        else {

            // note that the DqmcIO class should be a friend class of any derived model class
            // to get access to the bosonic fields member

            // ---------------------------------  Repulsive Hubbard model  ------------------------------------
            if ( const auto repulsive_hubbard = dynamic_cast<const Model::RepulsiveHubbard*>(&model);
                repulsive_hubbard != nullptr ) {
                // output current configuration of auxiliary bosonic fields
                // for repulsive hubbard model, they are ising-like.
                boost::format fmt_fields_info("%| 20d|%| 20d|");
                boost::format fmt_fields("%| 20d|%| 20d|%| 20.1f|");
                const int time_size = repulsive_hubbard->m_bosonic_field.rows();
                const int space_size = repulsive_hubbard->m_bosonic_field.cols();
                
                ostream << fmt_fields_info % time_size % space_size << std::endl;
                for ( auto t = 0; t < time_size; ++t ) {
                    for ( auto i = 0; i < space_size; ++i ) {
                        ostream << fmt_fields % t % i % repulsive_hubbard->m_bosonic_field(t,i) << std::endl;
                    }
                }
            }
            
            // ---------------------------------  Attractive Hubbard model  -----------------------------------
            else if ( const auto attractive_hubbard = dynamic_cast<const Model::AttractiveHubbard*>(&model); 
                attractive_hubbard != nullptr ) {
                // output current configuration of auxiliary bosonic fields
                // for attractive hubbard model, they are ising-like.
                boost::format fmt_fields_info("%| 20d|%| 20d|");
                boost::format fmt_fields("%| 20d|%| 20d|%| 20.1f|");
                const int time_size = attractive_hubbard->m_bosonic_field.rows();
                const int space_size = attractive_hubbard->m_bosonic_field.cols();
                
                ostream << fmt_fields_info % time_size % space_size << std::endl;
                for ( auto t = 0; t < time_size; ++t ) {
                    for ( auto i = 0; i < space_size; ++i ) {
                        ostream << fmt_fields % t % i % attractive_hubbard->m_bosonic_field(t,i) << std::endl;
                    }
                }
            }

            // other model types, raising errors
            else {
                std::cerr << "QuantumMonteCarlo::DqmcIO::output_bosonic_fields(): "
                          << "undefined model type." << std::endl;
                exit(1);
            }
        }
    }


    void DqmcIO::read_bosonic_fields_from_file ( const std::string& filename, ModelBase& model )
    {   
        std::ifstream infile(filename, std::ios::in);

        // check whether the ifstream works well
        if ( !infile.is_open() ) {
            std::cerr << "QuantumMonteCarlo::DqmcIO::read_bosonic_fields_from_file(): "
                      << "fail to open file \'" << filename << "\'." << std::endl;
            exit(1);
        }

        // temporary parameters
        std::string line;
        std::vector<std::string> data;

        // note that the DqmcIO class should be a friend class of any derived model class
        // to get access to the bosonic fields member

        // ---------------------------------  Repulsive Hubbard model  ------------------------------------
        if ( auto repulsive_hubbard = dynamic_cast<Model::RepulsiveHubbard*>(&model);
            repulsive_hubbard != nullptr ) {

            // consistency check of the model parameters
            // read the first line which containing the model information
            getline(infile, line);
            boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
            data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));

            const int time_size = boost::lexical_cast<int>(data[0]);
            const int space_size = boost::lexical_cast<int>(data[1]); 
            if (   ( time_size != repulsive_hubbard->m_bosonic_field.rows() ) 
                || ( space_size != repulsive_hubbard->m_bosonic_field.cols() ) ) {
                std::cerr << "QuantumMonteCarlo::DqmcIO::read_bosonic_fields_from_file(): "
                          << "inconsistency between model settings and input configs (time or space size). " 
                          << std::endl;
                exit(1);
            }

            // read in the configurations of auxiliary fields
            int time_point, space_point;
            while( getline(infile, line) ) {
                boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
                data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
                time_point = boost::lexical_cast<int>(data[0]);
                space_point = boost::lexical_cast<int>(data[1]);
                repulsive_hubbard->m_bosonic_field(time_point, space_point) = boost::lexical_cast<double>(data[2]);
            }
            // close the file stream
            infile.close();
        }

        // ---------------------------------  Attractive Hubbard model  -----------------------------------
        else if ( auto attractive_hubbard = dynamic_cast<Model::AttractiveHubbard*>(&model);
            attractive_hubbard != nullptr ) {

            // consistency check of the model parameters
            // read the first line which containing the model information
            getline(infile, line);
            boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
            data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));

            const int time_size = boost::lexical_cast<int>(data[0]);
            const int space_size = boost::lexical_cast<int>(data[1]); 
            if (   ( time_size != attractive_hubbard->m_bosonic_field.rows() ) 
                || ( space_size != attractive_hubbard->m_bosonic_field.cols() ) ) {
                std::cerr << "QuantumMonteCarlo::DqmcIO::read_bosonic_fields_from_file(): "
                          << "inconsistency between model settings and input configs (time or space size). " 
                          << std::endl;
                exit(1);
            }

            // read in the configurations of auxiliary fields
            int time_point, space_point;
            while( getline(infile, line) ) {
                boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
                data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
                time_point = boost::lexical_cast<int>(data[0]);
                space_point = boost::lexical_cast<int>(data[1]);
                attractive_hubbard->m_bosonic_field(time_point, space_point) = boost::lexical_cast<double>(data[2]);
            }
            // close the file stream
            infile.close();
        }

        // other model types, raising errors
        else {
            // close the file stream
            infile.close();
            std::cerr << "QuantumMonteCarlo::DqmcIO::read_bosonic_fields_from_file(): "
                      << "undefined model type." << std::endl;
            exit(1);
        }
    }



} // namespace QuantumMonteCarlo

#endif // DQMC_IO_H
