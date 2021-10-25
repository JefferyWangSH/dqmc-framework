#include "measure_data.h"
#include "eigen_boost_serialization.hpp"


template<class Archive>
void Measure::MeasureData::serialize(Archive & ar, const unsigned int version) {
    // TODO
};

Measure::MeasureData::MeasureData(const int &size_of_bin) {
    this->set_size_of_bin(size_of_bin);
}

int Measure::MeasureData::counts() const {
    return this->_count;
}

int Measure::MeasureData::size_of_bin() const {
    assert( this->_size_of_bin == this->_bins.size() );
    return this->_size_of_bin;
}

double Measure::MeasureData::mean_value() const {
    return this->_mean_value;
}

double Measure::MeasureData::error_bar() const {
    return this->_error_bar;
}

double Measure::MeasureData::tmp_value() const {
    return this->_tmp_data;
}

double& Measure::MeasureData::tmp_value() {
    return this->_tmp_data;
}

Eigen::VectorXd& Measure::MeasureData::bin_data() {
    return this->_bins;
}

const Eigen::VectorXd &Measure::MeasureData::bin_data() const {
    return this->_bins;
}

void Measure::MeasureData::set_size_of_bin(const int &size_of_bin) {
    this->_size_of_bin = size_of_bin;
    this->_bins.resize(size_of_bin);
}

void Measure::MeasureData::clear() {
    this->_mean_value = 0.0;
    this->_error_bar = 0.0;
}

void Measure::MeasureData::clear_temporary() {
    this->_tmp_data = 0.0;
    this->_count = 0;
}

void Measure::MeasureData::clear_bin_data() {
    this->_bins.setZero();
}

void Measure::MeasureData::calculate_mean_value() {
    this->_mean_value = this->_bins.sum() / this->size_of_bin();
}

void Measure::MeasureData::calculate_error_bar() {
    this->_error_bar = this->_bins.array().square().sum() / this->size_of_bin();
    this->_error_bar = pow(this->_error_bar - pow(this->_mean_value, 2), 0.5) / pow(this->size_of_bin()-1, 0.5);
}

void Measure::MeasureData::analyse() {
    this->calculate_mean_value();
    this->calculate_error_bar();
}
