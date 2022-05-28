#ifndef OBSERABLE_H
#define OBSERABLE_H
#pragma once

/**
  *  This head file includes the abstract base class Observable::ObservableBase 
  *  and its template derived class Observable::Observable<ObsType> class, which are
  *  designed for the measurements of physical observables in dqmc simualtions.
  *  Support observable types of scalar, vector and matrix kinds.
  */

#include <vector>
#include <string>
#include <functional>
#include <numeric>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


// forward declaration
namespace QuantumMonteCarlo { class DqmcWalker; }
namespace Model { class ModelBase; }
namespace Lattice { class LatticeBase; }
namespace Measure { class MeasureHandler; }


namespace Observable {

    // data types of observables
    using ScalarType = double;
    using VectorType = Eigen::VectorXd;
    using MatrixType = Eigen::MatrixXd;

    
    // --------------------------- Abstract base class Observable::ObservableBase ----------------------------
    // this class should not be instantiated in any case 
    // it only serves as a pointer to its derived class
    class ObservableBase {
        protected:
            ObservableBase() = default;
            
        public:
            virtual ~ObservableBase(){};
    };

    
    // -----------------------  Derived template class Observable::Observable<ObsType> -----------------------
    template<typename ObsType> class Observable : public ObservableBase {
        private:
            
            using DqmcWalker = QuantumMonteCarlo::DqmcWalker;
            using ModelBase = Model::ModelBase;
            using LatticeBase = Lattice::LatticeBase;
            using MeasureHandler = Measure::MeasureHandler;
            using ObsMethod = void( Observable<ObsType>&,
                                    const MeasureHandler&, 
                                    const DqmcWalker&, 
                                    const ModelBase&,
                                    const LatticeBase& );

            ObsType m_mean_value{};
            ObsType m_error_bar{};
            ObsType m_tmp_value{};
            ObsType m_zero_elem{};

            std::string m_name{};
            int m_count{0};
            int m_bin_num{0};
            std::vector<ObsType> m_bin_data{};

            // user-defined method of measurements
            std::function<ObsMethod> m_method{};

        
        public:
           
            Observable() = default;
            
            explicit Observable(int bin_num) { this->set_number_of_bins(bin_num); }

            // overload operator ++
            int operator++() { return ++this->m_count; }


            // --------------------------------- Interface functions -----------------------------------
            
            int& counts() { return this->m_count; }
            int  counts() const { return this->m_count; }
            int  bin_num() const { return this->m_bin_num; }
            std::string name() const { return this->m_name; }
            
            const ObsType& zero_element() const { return this->m_zero_elem; } 
            const ObsType& mean_value() const { return this->m_mean_value; }
            const ObsType& error_bar() const { return this->m_error_bar; }
            
            const ObsType& tmp_value() const { return this->m_tmp_value; }
            ObsType& tmp_value() { return this->m_tmp_value; }

            const ObsType& bin_data(int bin) const {
                assert( bin >= 0 && bin < this->m_bin_num );
                return this->m_bin_data[bin];
            }

            ObsType& bin_data(int bin) {
                assert( bin >= 0 && bin < this->m_bin_num );
                return this->m_bin_data[bin];
            }


            // ----------------------------- Set up parameters and methods -----------------------------
            
            void set_number_of_bins(const int& bin_num) { this->m_bin_num = bin_num; }
            void set_zero_element(const ObsType& zero_elem) { this->m_zero_elem = zero_elem; }
            void set_observable_name(const std::string& name) { this->m_name = name; }
            void add_method(const std::function<ObsMethod>& method) { this->m_method = method; }


            // ------------------------------- Other member functions ----------------------------------
            
            // perform one step of measurment
            void measure( const MeasureHandler& meas_handler,
                          const DqmcWalker& walker, 
                          const ModelBase& model, 
                          const LatticeBase& lattice )
            { 
                this->m_method( *this, meas_handler, walker, model, lattice ); 
            }

            // allocate memory
            void allocate() {
                this->m_mean_value = this->m_zero_elem;
                this->m_error_bar = this->m_zero_elem;
                this->m_tmp_value = this->m_zero_elem;

                std::vector<ObsType>().swap(this->m_bin_data);
                this->m_bin_data.reserve(this->m_bin_num);
                for (int i = 0; i < this->m_bin_num; ++i) {
                    this->m_bin_data.emplace_back(this->m_zero_elem);
                }
            }

            void clear_stats() {
                this->m_mean_value = this->m_zero_elem;
                this->m_error_bar = this->m_zero_elem;
            }

            void clear_temporary() {
                this->m_tmp_value = this->m_zero_elem;
                this->m_count = 0;
            }

            void clear_bin_data() {
                for (auto& bin_data : this->m_bin_data) {
                    bin_data = this->m_zero_elem;
                }
            }

            // perform the analysis
            void analyse() {
                this->clear_stats();
                this->calculate_mean_value();
                this->calculate_error_bar();
            }


        private:

            // calculating mean value of measurements
            void calculate_mean_value() {
                this->m_mean_value = std::accumulate(this->m_bin_data.begin(), this->m_bin_data.end(), this->m_zero_elem);
                this->m_mean_value /= this->bin_num();
            }
            
            // the error bar is calculated according to the specific data type of observables,
            // which are defined in cpp file as specialized template member functions.
            void calculate_error_bar();
    };


    // declaration of specialized template member functions
    template<> void Observable<ScalarType>::calculate_error_bar();
    template<> void Observable<VectorType>::calculate_error_bar();
    template<> void Observable<MatrixType>::calculate_error_bar();


    // some aliases
    using ScalarObs = Observable<ScalarType>;
    using VectorObs = Observable<VectorType>;
    using MatrixObs = Observable<MatrixType>;


} // namespace Observable

#endif // OBSERVABLE_H
