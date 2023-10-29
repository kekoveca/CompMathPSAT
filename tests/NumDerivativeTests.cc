#include <gtest/gtest.h>
#include "../src/NumDerivative.hpp"
#include <iostream>
#include <fstream>

const std::string dataPath = "/Users/kekoveca/Documents/MIPT/CompMath/post_proccessing/NumDerivativeFolder";

TEST(NumDerivativeTests, simpleBasicFirst){
    std::array<double, 2> points = {-1, 1};
    auto coeffs = calcDerivativeCoef<double, 2, 1>(points);
    ASSERT_NEAR(coeffs.centralCoef, 0., 1e-15);
    ASSERT_NEAR(coeffs.otherCoefs[0], -0.5, 1e-15);
    ASSERT_NEAR(coeffs.otherCoefs[1], 0.5, 1e-15);
}

TEST(NumDerivativeTests, simpleBasicSecond){
    std::array<double, 2> points = {1, 2};
    auto coeffs = calcDerivativeCoef<double, 2, 1>(points);
    ASSERT_NEAR(coeffs.otherCoefs[0], 2, 1e-15);
    ASSERT_NEAR(coeffs.otherCoefs[1], -0.5, 1e-15);
}

TEST(NumDerivativeTests, mainBasic){
    std::ofstream fileOut;
    fileOut.open(dataPath + "/errorsBasic.txt");

    constexpr double x = 1.;

    constexpr unsigned int L = 1;

    const double exp0 = std::exp(x); 

    constexpr unsigned int N3 = 3;
    constexpr unsigned int N4 = 4;
    constexpr unsigned int N5 = 5;
    
    std::array<double, N3> pointsN3 = {-1., 1., 2.};
    std::array<double, N4> pointsN4 = {-2., -1., 1., 2.};
    std::array<double, N5> pointsN5 = {-3., -2., -1., 1., 2.};

    std::array<double, 16> h;

    h[0] = 1;

    for(unsigned int i = 0; i < h.size() - 1; i++) {
        h[i + 1] = h[i] * 0.1;
    }

    auto coeffsN3 = calcDerivativeCoef<double, N3, L>(pointsN3);

    for(auto&& iter : h) {
        double DN3 = coeffsN3.centralCoef * exp0;
        for(unsigned int i = 0; i < pointsN3.size(); i++){
            DN3 += coeffsN3.otherCoefs[i] * std::exp(x + pointsN3[i] * iter);
        }
        fileOut << std::abs(DN3 / iter - exp0)  <<  ",";
    }
    fileOut << std::endl;

    auto coeffsN4 = calcDerivativeCoef<double, N4, L>(pointsN4);

    for(auto&& iter : h) {
        double DN4 = coeffsN4.centralCoef * exp0;
        for(unsigned int i = 0; i < pointsN4.size(); i++){
            DN4 += coeffsN4.otherCoefs[i] * std::exp(x + pointsN4[i] * iter);
        }
        fileOut << std::abs(DN4 / iter - exp0)  <<  ",";
    }
    fileOut << std::endl;

    auto coeffsN5 = calcDerivativeCoef<double, N5, L>(pointsN5);

    for(auto&& iter : h) {
        double DN5 = coeffsN5.centralCoef * exp0;
        for(unsigned int i = 0; i < pointsN5.size(); i++){
            DN5 += coeffsN5.otherCoefs[i] * std::exp(x + pointsN5[i] * iter);
        }
        fileOut << std::abs(DN5 / iter - exp0)  <<  ",";
    }
    fileOut << std::endl;
    fileOut.close();
}

TEST(NumDerivativeTests, simpleAdvanced){
    std::array<double, 2> points = {-1, 1};
    auto coeffs = calcDerivativeCoef<double, 2, 2>(points);
    ASSERT_NEAR(coeffs.centralCoef, -2., 1e-15);
    ASSERT_NEAR(coeffs.otherCoefs[0], 1., 1e-15);
    ASSERT_NEAR(coeffs.otherCoefs[1], 1., 1e-15);
}

TEST(NumDerivativeTests, mainAdvanced){
    std::ofstream fileOut;
    fileOut.open(dataPath + "/errorsAdvanced.txt");

    constexpr double x = 1.;

    constexpr unsigned int L = 2;

    const double exp0 = std::exp(x); 

    constexpr unsigned int N3 = 3;
    constexpr unsigned int N4 = 4;
    constexpr unsigned int N5 = 5;
    
    std::array<double, N3> pointsN3 = {-1., 1., 2.};
    std::array<double, N4> pointsN4 = {-1., 1., 2., 3.};
    std::array<double, N5> pointsN5 = {-2., -1., 1., 2., 3.};

    std::array<double, 16> h;

    h[0] = 1;

    for(unsigned int i = 0; i < h.size() - 1; i++) {
        h[i + 1] = h[i] * 0.1;
    }

    auto coeffsN3 = calcDerivativeCoef<double, N3, L>(pointsN3);

    for(auto&& iter : h) {
        double DN3 = coeffsN3.centralCoef * exp0;
        for(unsigned int i = 0; i < pointsN3.size(); i++){
            DN3 += coeffsN3.otherCoefs[i] * std::exp(x + pointsN3[i] * iter);
        }
        fileOut << std::abs(DN3 / std::pow(iter, L) - exp0)  <<  ",";
    }
    fileOut << std::endl;

    auto coeffsN4 = calcDerivativeCoef<double, N4, L>(pointsN4);

    for(auto&& iter : h) {
        double DN4 = coeffsN4.centralCoef * exp0;
        for(unsigned int i = 0; i < pointsN4.size(); i++){
            DN4 += coeffsN4.otherCoefs[i] * std::exp(x + pointsN4[i] * iter);
        }
        fileOut << std::abs(DN4 / std::pow(iter, L) - exp0)  <<  ",";
    }
    fileOut << std::endl;

    auto coeffsN5 = calcDerivativeCoef<double, N5, L>(pointsN5);

    for(auto&& iter : h) {
        double DN5 = coeffsN5.centralCoef * exp0;
        for(unsigned int i = 0; i < pointsN5.size(); i++){
            DN5 += coeffsN5.otherCoefs[i] * std::exp(x + pointsN5[i] * iter);
        }
        fileOut << std::abs(DN5 / std::pow(iter, L) - exp0)  <<  ",";
    }
    fileOut << std::endl;
    fileOut.close();
}
