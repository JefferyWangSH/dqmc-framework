#include "measure_data.h"

/* implementations of specialized template functions */

template<> void Measure::MeasureData<double>::calculate_error_bar() {
    for (auto bin_data : this->_bins) {
        this->_error_bar += std::pow(bin_data, 2);
    }
    this->_error_bar /= this->size_of_bin();
    this->_error_bar = std::sqrt(this->_error_bar - std::pow(this->_mean_value,2)) / std::sqrt(this->size_of_bin()-1);
}

template<> void Measure::MeasureData<Eigen::VectorXd>::calculate_error_bar() {
    for (auto bin_data : this->_bins) {
        this->_error_bar += bin_data.array().square().matrix();
    }
    this->_error_bar /= this->size_of_bin();
    this->_error_bar = (this->_error_bar.array() - this->_mean_value.array().square()).sqrt().matrix() / std::sqrt(this->size_of_bin()-1);
}

template<> void Measure::MeasureData<Eigen::MatrixXd>::calculate_error_bar() {
    for (auto bin_data : this->_bins) {
        this->_error_bar += bin_data.array().square().matrix();
    }
    this->_error_bar /= this->size_of_bin();
    this->_error_bar = (this->_error_bar.array() - this->_mean_value.array().square()).sqrt().matrix() / std::sqrt(this->size_of_bin()-1);
}
