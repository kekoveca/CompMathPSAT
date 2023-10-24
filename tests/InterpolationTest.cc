#include <gtest/gtest.h>
#include "../src/Interpolation.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>

const std::string dataPath = "/Users/kekoveca/Documents/MIPT/CompMath/post_proccessing/InterpolationFolder";

template<typename xType, unsigned int N>
std::array<xType, N> chebyshevZeros(xType a, xType b){
    std::array<xType, N> zeros;
    for(unsigned int i = 0; i < N; i++){
        zeros[i] = (b + a) / xType(2) + (b - a) * std::cos(M_PI * (2 * i + 1)/(2 * N)) / xType(2);
    }
    return zeros;
}


TEST(InterpolantTests, simpleEigenTest){
    const Eigen::Array<double, 3, 1> points = {0., 1, 2};
    const Eigen::Array<double, 3, 1> values = {1, 3, 9};

    EigenNewtonInterpolator<double, double, 3> interpolator(points, values);

    auto diffs = interpolator.getsplitDifferencies();

    ASSERT_EQ(diffs[0], 1.);
    ASSERT_EQ(diffs[1], 2.);
    ASSERT_EQ(diffs[2], 2.);

    ASSERT_EQ(interpolator.interpolate(1), 3);
}

TEST(InterpolantTests, simpleTest){
    const std::array<double, 3> points = {0., 1, 2};
    const std::array<double, 3> values = {1, 3, 9};

    NewtonInterpolator<double, double, 3> interpolator(points, values);

    auto diffs = interpolator.getsplitDifferencies();

    ASSERT_EQ(diffs[0], 1.);
    ASSERT_EQ(diffs[1], 2.);
    ASSERT_EQ(diffs[2], 2.);
}


TEST(InterpolantTests, Test){
    std::ofstream fileOutN3;
    std::ofstream fileOutN4;
    std::ofstream fileOutN5;

    fileOutN3.open(dataPath + "/errN3.txt");
    fileOutN4.open(dataPath + "/errN4.txt");
    fileOutN5.open(dataPath + "/errN5.txt");

    constexpr int power = 6;

    // N = 3

    Eigen::ArrayXXd pointsMatrix3(3, 6);
    Eigen::ArrayXXd valuesMatrix3(3, 6);
    Eigen::ArrayXXd pointsMatrixErr3(1000, 6);
    Eigen::ArrayXXd valuesMatrixErr3(1000, 6);

    for(int i = 0; i < power; i++) {
        pointsMatrix3.col(i) = Eigen::ArrayXd::LinSpaced(3, 0., std::pow(2., 1 - i));
        valuesMatrix3.col(i) = pointsMatrix3.col(i).exp();
    }

    for(int i = 0; i < power; i++) {
        pointsMatrixErr3.col(i) = Eigen::ArrayXd::LinSpaced(1000, 0., std::pow(2., 1 - i));
        valuesMatrixErr3.col(i) = pointsMatrixErr3.col(i).exp();
    }

    for(int i = 0; i < power; i++) {
        double err = 0;
        EigenNewtonInterpolator<double, double, 3>interpolator(pointsMatrix3.col(i), valuesMatrix3.col(i));
        for(int j = 0; j < 1000; j++){
            auto fInterpolated = interpolator.interpolate(pointsMatrixErr3.col(i)(j));
            err = std::max(std::abs(valuesMatrixErr3.col(i)(j) - fInterpolated), err);
        }
        fileOutN3 << err << "\n";
    }

    // N = 4

    Eigen::ArrayXXd pointsMatrix4(4, 6);
    Eigen::ArrayXXd valuesMatrix4(4, 6);
    Eigen::ArrayXXd pointsMatrixErr4(1000, 6);
    Eigen::ArrayXXd valuesMatrixErr4(1000, 6);

    for(int i = 0; i < power; i++) {
        pointsMatrix4.col(i) = Eigen::ArrayXd::LinSpaced(4, 0., std::pow(2., 1 - i));
        valuesMatrix4.col(i) = pointsMatrix4.col(i).exp();
    }

    for(int i = 0; i < power; i++) {
        pointsMatrixErr4.col(i) = Eigen::ArrayXd::LinSpaced(1000, 0., std::pow(2., 1 - i));
        valuesMatrixErr4.col(i) = pointsMatrixErr4.col(i).exp();
    }

    for(int i = 0; i < power; i++) {
        double err = 0;
        EigenNewtonInterpolator<double, double, 4>interpolator(pointsMatrix4.col(i), valuesMatrix4.col(i));
        for(int j = 0; j < 1000; j++){
            auto fInterpolated = interpolator.interpolate(pointsMatrixErr4.col(i)(j));
            err = std::max(std::abs(valuesMatrixErr4.col(i)(j) - fInterpolated), err);
        }
        fileOutN4 << err << "\n";
    }

    // N = 5

    Eigen::ArrayXXd pointsMatrix5(5, 6);
    Eigen::ArrayXXd valuesMatrix5(5, 6);
    Eigen::ArrayXXd pointsMatrixErr5(1000, 6);
    Eigen::ArrayXXd valuesMatrixErr5(1000, 6);

    for(int i = 0; i < power; i++) {
        pointsMatrix5.col(i) = Eigen::ArrayXd::LinSpaced(5, 0., std::pow(2., 1 - i));
        valuesMatrix5.col(i) = pointsMatrix5.col(i).exp();
    }

    for(int i = 0; i < power; i++) {
        pointsMatrixErr5.col(i) = Eigen::ArrayXd::LinSpaced(1000, 0., std::pow(2., 1 - i));
        valuesMatrixErr5.col(i) = pointsMatrixErr5.col(i).exp();
    }

    for(int i = 0; i < power; i++) {
        double err = 0;
        EigenNewtonInterpolator<double, double, 5>interpolator(pointsMatrix5.col(i), valuesMatrix5.col(i));
        for(int j = 0; j < 1000; j++){
            auto fInterpolated = interpolator.interpolate(pointsMatrixErr5.col(i)(j));
            err = std::max(std::abs(valuesMatrixErr5.col(i)(j) - fInterpolated), err);
        }
        fileOutN5 << err << "\n";
    }
}