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
namespace Measure { class MeasureHandler; }
namespace Model { class ModelBase; }
namespace Lattice { class LatticeBase; }


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
            
            using MeasureHandler = Measure::MeasureHandler;
            using ModelBase = Model::ModelBase;
            using LatticeBase = Lattice::LatticeBase;
            using ObsMethod = void( Observable<ObsType>&, 
                                    const MeasureHandler&, 
                                    const ModelBase&,
                                    const LatticeBase& );

            ObsType m_mean_value{};
            ObsType m_error_bar{};
            ObsType m_tmp_data{};
            ObsType m_zero_elem{};

            std::string m_name{};
            int m_count{0};
            int m_size_of_bin{0};
            std::vector<ObsType> m_bin_data{};

            // user-defined method of measurements
            std::function<ObsMethod> m_method{};

        
        public:
           
            Observable() = default;
            
            explicit Observable(int size_of_bin) { this->set_size_of_bin(size_of_bin); }

            // overload operator ++
            int operator++() { return ++this->m_count; }


            // --------------------------------- Interface functions -----------------------------------
            
            int counts() const { return this->m_count; }
            int size_of_bin() const { return this->m_size_of_bin; }
            std::string name() const { return this->m_name; }
            
            const ObsType& zero_element() const { return this->m_zero_elem; } 
            const ObsType& mean_value() const { return this->m_mean_value; }
            const ObsType& error_bar() const { return this->m_error_bar; }
            
            const ObsType& tmp_value() const { return this->tmp_value; }
            ObsType& tmp_value() { return this->tmp_value; }

            const std::vector<ObsType>& bin_data() const { return this->m_bin_data; }
            std::vector<ObsType>& bin_data() { return this->m_bin_data; }


            // ----------------------------- Set up parameters and methods -----------------------------
            
            void set_size_of_bin(const int& size_of_bin) { this->m_size_of_bin = size_of_bin; }
            void set_zero_element(const ObsType& zero_elem) { this->m_zero_elem = zero_elem; }
            void set_observable_name(const std::string& name) { this->m_name = name; }
            void add_method(const std::function<ObsMethod>& method) { this->m_method = method; }


            // ------------------------------- Other member functions ----------------------------------
            
            // perform one step of measurment
            void measure(   const MeasureHandler& meas_handler, 
                            const ModelBase& model, 
                            const LatticeBase& lattice  )
            { 
                this->m_method(*this, meas_handler, model, lattice); 
            }

            // allocate memory
            void allocate() {
                this->m_mean_value = this->m_zero_elem;
                this->m_error_bar = this->m_zero_elem;
                this->m_tmp_data = this->m_zero_elem;

                std::vector<ObsType>().swap(this->m_bin_data);
                this->m_bin_data.reserve(this->m_size_of_bin);
                for (int i = 0; i < this->m_size_of_bin; ++i) {
                    this->m_bin_data.emplace_back(this->m_zero_elem);
                }
            }

            void clear_stats() {
                this->m_mean_value = this->m_zero_elem;
                this->m_error_bar = this->m_zero_elem;
            }

            void clear_temporary() {
                this->m_tmp_data = this->m_zero_elem;
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
                this->m_mean_value /= this->size_of_bin();
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
