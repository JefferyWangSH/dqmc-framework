#ifndef DQMC_HUBBARD_MEASUREDATA_H
#define DQMC_HUBBARD_MEASUREDATA_H
#pragma once

/**
  *  This head file includes data structure `MeasureData`,
  *  designed for the measurements of physical observables.
  *  Support data types of double, Eigen::VectorXd and Eigen::MatrixXd using templates.
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include <vector>


namespace Measure {

    template<class DataStructure>
    class MeasureData {
    private:
        using BinStructure = std::vector<DataStructure>;
        DataStructure _mean_value{};
        DataStructure _error_bar{};
        DataStructure _tmp_data{};
        DataStructure _zero_elem{};

        int _count{0};
        int _size_of_bin{0};
        BinStructure _bins{};

    public:
        MeasureData() = default;
        explicit MeasureData(const int& size_of_bin);

        int operator++() { return ++this->_count; }

        int counts() const;

        int size_of_bin() const;

        const DataStructure& zero_element() const; 

        const DataStructure& mean_value() const;

        const DataStructure& error_bar() const;

        const DataStructure& tmp_value() const;

        DataStructure& tmp_value();

        const BinStructure& bin_data() const;

        BinStructure& bin_data();

        void set_size_of_bin(const int &size_of_bin);

        void set_zero_element(const DataStructure& zero_elem);

        void allocate();

        void clear();

        void clear_temporary();

        void clear_bin_data();

        void calculate_mean_value();

        void calculate_error_bar();

        void analyse();
    };

    // declaration of specialized template functions
    template<> void Measure::MeasureData<double>::calculate_error_bar();
    template<> void Measure::MeasureData<Eigen::VectorXd>::calculate_error_bar();
    template<> void Measure::MeasureData<Eigen::MatrixXd>::calculate_error_bar();



    /* Implementations of member functions with templates */

    template<class DataStructure>
    Measure::MeasureData<DataStructure>::MeasureData(const int &size_of_bin) {
        this->set_size_of_bin(size_of_bin);
    }

    template<class DataStructure>
    int Measure::MeasureData<DataStructure>::counts() const {
        return this->_count;
    }

    template<class DataStructure>
    int Measure::MeasureData<DataStructure>::size_of_bin() const {
        assert( this->_size_of_bin == this->_bins.size() );
        return this->_size_of_bin;
    }

    template<class DataStructure>
    const DataStructure& Measure::MeasureData<DataStructure>::zero_element() const {
        return this->_zero_elem;
    }

    template<class DataStructure>
    const DataStructure& Measure::MeasureData<DataStructure>::mean_value() const {
        return this->_mean_value;
    }

    template<class DataStructure>
    const DataStructure& Measure::MeasureData<DataStructure>::error_bar() const {
        return this->_error_bar;
    }

    template<class DataStructure>
    const DataStructure& Measure::MeasureData<DataStructure>::tmp_value() const {
        return this->_tmp_data;
    }

    template<class DataStructure>
    DataStructure& Measure::MeasureData<DataStructure>::tmp_value() {
        return this->_tmp_data;
    }

    template<class DataStructure>
    std::vector<DataStructure>& Measure::MeasureData<DataStructure>::bin_data() {
        return this->_bins;
    }

    template<class DataStructure>
    const std::vector<DataStructure>& Measure::MeasureData<DataStructure>::bin_data() const {
        return this->_bins;
    }

    template<class DataStructure>
    void Measure::MeasureData<DataStructure>::set_zero_element(const DataStructure &zero_elem) {
        this->_zero_elem = zero_elem;
    }

    template<class DataStructure>
    void Measure::MeasureData<DataStructure>::set_size_of_bin(const int &size_of_bin) {
        this->_size_of_bin = size_of_bin;
    }

    template<class DataStructure>
    void Measure::MeasureData<DataStructure>::allocate() {
        this->_mean_value = this->_zero_elem;
        this->_error_bar = this->_zero_elem;
        this->_tmp_data = this->_zero_elem;

        this->_bins.clear();
        this->_bins.reserve(this->_size_of_bin);
        for (int i = 0; i < this->_size_of_bin; ++i) {
            this->_bins.emplace_back(this->_zero_elem);
        }
    }

    template<class DataStructure>
    void Measure::MeasureData<DataStructure>::clear() {
        this->_mean_value = this->_zero_elem;
        this->_error_bar = this->_zero_elem;
    }

    template<class DataStructure>
    void Measure::MeasureData<DataStructure>::clear_temporary() {
        this->_tmp_data = this->_zero_elem;
        this->_count = 0;
    }

    template<class DataStructure>
    void Measure::MeasureData<DataStructure>::clear_bin_data() {
        for (auto bin_data : this->_bins) {
            bin_data = this->_zero_elem;
        }
    }

    template<class DataStructure>
    void Measure::MeasureData<DataStructure>::calculate_mean_value() {
        this->_mean_value = std::accumulate(this->_bins.begin(), this->_bins.end(), this->_zero_elem) / this->size_of_bin();
    }

    template<class DataStructure>
    void Measure::MeasureData<DataStructure>::analyse() {
        this->clear();
        this->calculate_mean_value();
        this->calculate_error_bar();
    }
}

#endif //DQMC_HUBBARD_MEASUREDATA_H
