#ifndef DQMC_HUBBARD_MEASUREDATA_H
#define DQMC_HUBBARD_MEASUREDATA_H
#pragma once

/**
  *  This is head file includes data structure `MeasureData` designed for the measurements of physical observables.
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

namespace Measure {

    class MeasureData {
    private:
        double _mean_value = 0.0;
        double _error_bar = 0.0;
        double _tmp_data = 0.0;

        int _size_of_bin = 0.0;
        Eigen::VectorXd _bins{};

    public:
        MeasureData() = default;
        explicit MeasureData(const int& size_of_bin);

        double mean_value() const;

        double error_bar() const;

        double tmp_value() const;

        double& tmp_value();

        int size_of_bin() const;

        const Eigen::VectorXd& bin_data() const;

        Eigen::VectorXd& bin_data();

        void set_size_of_bin(const int &size_of_bin);

        void clear();

        void clear_temporary();

        void clear_bin_data();

        void calculate_mean_value();

        void calculate_error_bar();

        void analyse();
    };
}

#endif //DQMC_HUBBARD_MEASUREDATA_H
