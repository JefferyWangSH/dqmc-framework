#include "measure/observable.h"


namespace Observable {

    // implementations of specialized template functions

    template<> void Observable<ScalarType>::calculate_error_bar() {
        for (const auto& bin_data : this->m_bin_data) {
            this->m_error_bar += std::pow(bin_data, 2);
        }
        this->m_error_bar /= this->size_of_bin();
        this->m_error_bar = std::sqrt(this->m_error_bar - std::pow(this->m_mean_value,2)) / std::sqrt(this->size_of_bin()-1);
    }

    template<> void Observable<VectorType>::calculate_error_bar() {
        for (const auto& bin_data : this->m_bin_data) {
            this->m_error_bar += bin_data.array().square().matrix();
        }
        this->m_error_bar /= this->size_of_bin();
        this->m_error_bar = (this->m_error_bar.array() - this->m_mean_value.array().square()).sqrt().matrix() / std::sqrt(this->size_of_bin()-1);
    }

    template<> void Observable<MatrixType>::calculate_error_bar() {
        for (const auto& bin_data : this->m_bin_data) {
            this->m_error_bar += bin_data.array().square().matrix();
        }
        this->m_error_bar /= this->size_of_bin();
        this->m_error_bar = (this->m_error_bar.array() - this->m_mean_value.array().square()).sqrt().matrix() / std::sqrt(this->size_of_bin()-1);
    }


} // namespace Observable