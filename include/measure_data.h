#ifndef DQMC_HUBBARD_MEASUREDATA_H
#define DQMC_HUBBARD_MEASUREDATA_H
#pragma once

/**
  *  This head file includes data structure `MeasureData` designed for the measurements of physical observables.
  *  and the seiraliztion of class are realized using boost.
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include <boost/serialization/access.hpp>


namespace Measure {

    class MeasureData {
    private:
        double _mean_value = 0.0;
        double _error_bar = 0.0;
        double _tmp_data = 0.0;

        int _count = 0;

        int _size_of_bin = 0.0;
        Eigen::VectorXd _bins{};

    private:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version);

    public:
        MeasureData() = default;
        explicit MeasureData(const int& size_of_bin);

        int operator++() { return ++this->_count; }

        int counts() const;

        int size_of_bin() const;

        double mean_value() const;

        double error_bar() const;

        double tmp_value() const;

        double& tmp_value();

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
