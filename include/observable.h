#ifndef DQMC_HUBBARD_OBSERABLE_H
#define DQMC_HUBBARD_OBSERABLE_H
#pragma once

/**
  *  This head file includes `Measure::Observable` class
  *  designed for the measurements of physical observables.
  *  Support data types of double, Eigen::VectorXd and Eigen::MatrixXd using templates.
  */

#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <numeric>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


// forward declaration
namespace Model { class Hubbard; }
namespace Measure { class Measure; }

namespace Measure {

    template<class DataType>
    class Observable {
    private:
        using BinDataType = std::vector<DataType>;
        DataType _mean_value{};
        DataType _error_bar{};
        DataType _tmp_data{};
        DataType _zero_elem{};

        std::string _name{};
        int _count{0};
        int _size_of_bin{0};
        BinDataType _bins{};

        // user-defined method of measurements
        using MethodTemplate = void(Observable<DataType>&, Measure&, const Model::Hubbard&);
        std::function<MethodTemplate> _method{};
        

    public:
        Observable() = default;
        explicit Observable(const int& size_of_bin);

        int operator++() { return ++this->_count; }

        int counts() const;

        int size_of_bin() const;
        
        const std::string& name() const;

        const DataType& zero_element() const; 

        const DataType& mean_value() const;

        const DataType& error_bar() const;

        const DataType& tmp_value() const;

        DataType& tmp_value();

        const BinDataType& bin_data() const;

        BinDataType& bin_data();

        void set_size_of_bin(const int &size_of_bin);

        void set_zero_element(const DataType &zero_elem);

        void set_observable_name(const std::string &name);
        
        void add_method(const std::function<MethodTemplate> &method);

        void measure(Measure &measure, const Model::Hubbard &hubbard);

        void allocate();

        void clear_stats();

        void clear_temporary();

        void clear_bin_data();

        void calculate_mean_value();

        void calculate_error_bar();

        void analyse();
    };

    // declaration of specialized template functions
    template<> void Observable<double>::calculate_error_bar();
    template<> void Observable<Eigen::VectorXd>::calculate_error_bar();
    template<> void Observable<Eigen::MatrixXd>::calculate_error_bar();



    /* Implementations of member functions with templates */

    template<class DataType>
    Observable<DataType>::Observable(const int &size_of_bin) {
        this->set_size_of_bin(size_of_bin);
    }

    template<class DataType>
    int Observable<DataType>::counts() const {
        return this->_count;
    }

    template<class DataType>
    int Observable<DataType>::size_of_bin() const {
        assert( this->_size_of_bin == this->_bins.size() );
        return this->_size_of_bin;
    }

    template<class DataType>
    const std::string& Observable<DataType>::name() const {
        return this->_name;
    }

    template<class DataType>
    const DataType& Observable<DataType>::zero_element() const {
        return this->_zero_elem;
    }

    template<class DataType>
    const DataType& Observable<DataType>::mean_value() const {
        return this->_mean_value;
    }

    template<class DataType>
    const DataType& Observable<DataType>::error_bar() const {
        return this->_error_bar;
    }

    template<class DataType>
    const DataType& Observable<DataType>::tmp_value() const {
        return this->_tmp_data;
    }

    template<class DataType>
    DataType& Observable<DataType>::tmp_value() {
        return this->_tmp_data;
    }

    template<class DataType>
    std::vector<DataType>& Observable<DataType>::bin_data() {
        return this->_bins;
    }

    template<class DataType>
    const std::vector<DataType>& Observable<DataType>::bin_data() const {
        return this->_bins;
    }

    template<class DataType>
    void Observable<DataType>::set_zero_element(const DataType &zero_elem) {
        this->_zero_elem = zero_elem;
    }

    template<class DataType>
    void Observable<DataType>::set_size_of_bin(const int &size_of_bin) {
        this->_size_of_bin = size_of_bin;
    }

    template<class DataType>
    void Observable<DataType>::set_observable_name(const std::string &name) {
        this->_name = name;
    }

    template<class DataType>
    void Observable<DataType>::add_method(const std::function<MethodTemplate> &method) {
        this->_method = method;
    }

    template<class DataType>
    void Observable<DataType>::measure(Measure &measure, const Model::Hubbard &hubbard) {
        this->_method(*this, measure, hubbard);
    }

    template<class DataType>
    void Observable<DataType>::allocate() {
        this->_mean_value = this->_zero_elem;
        this->_error_bar = this->_zero_elem;
        this->_tmp_data = this->_zero_elem;

        std::vector<DataType>().swap(this->_bins);
        this->_bins.reserve(this->_size_of_bin);
        for (int i = 0; i < this->_size_of_bin; ++i) {
            this->_bins.emplace_back(this->_zero_elem);
        }
    }

    template<class DataType>
    void Observable<DataType>::clear_stats() {
        this->_mean_value = this->_zero_elem;
        this->_error_bar = this->_zero_elem;
    }

    template<class DataType>
    void Observable<DataType>::clear_temporary() {
        this->_tmp_data = this->_zero_elem;
        this->_count = 0;
    }

    template<class DataType>
    void Observable<DataType>::clear_bin_data() {
        for (auto bin_data : this->_bins) {
            bin_data = this->_zero_elem;
        }
    }

    template<class DataType>
    void Observable<DataType>::calculate_mean_value() {
        this->_mean_value = std::accumulate(this->_bins.begin(), this->_bins.end(), this->_zero_elem) / this->size_of_bin();
    }

    template<class DataType>
    void Observable<DataType>::analyse() {
        this->clear_stats();
        this->calculate_mean_value();
        this->calculate_error_bar();
    }

} // namespace Measure

#endif //DQMC_HUBBARD_OBSERVABLE_H
