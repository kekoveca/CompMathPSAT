#include <gtest/gtest.h>
#include "../src/NumDerivative.hpp"
#include <iostream>
#include <fstream>

const std::string dataPath = "/Users/kekoveca/Documents/MIPT/CompMath/post_proccessing/NumDerivativeFolder";

TEST(NumDerivativeTests, simple){
    std::array<double, 2> points = {-1, 1};
    auto coeffs = calcDerivativeCoef(points);
    ASSERT_NEAR(coeffs[0], 0, 1e-15);
    ASSERT_NEAR(coeffs[1], -0.5, 1e-15);
    ASSERT_NEAR(coeffs[2], 0.5, 1e-15);
}

TEST(NumDerivativeTests, main){
    std::ofstream fileOut;
    fileOut.open(dataPath + "/errors.txt");

    constexpr double x = 1.;

    const double exp0 = std::exp(x); 

    constexpr unsigned int N3 = 3;
    constexpr unsigned int N4 = 4;
    constexpr unsigned int N5 = 5;
    
    std::array<double, N3> pointsN3 = {1., 2., 3.};
    std::array<double, N4> pointsN4 = {1., 2., 3., 4.};
    std::array<double, N5> pointsN5 = {1., 2., 3., 4., 5.};

    std::array<double, 16> h;

    h[0] = 1;

    for(unsigned int i = 0; i < h.size() - 1; i++) {
        h[i + 1] = h[i] * 0.5;
    }

    auto coeffsN3 = calcDerivativeCoef(pointsN3);

    for(auto&& iter : h) {
        double DN3 = coeffsN3[0] * exp0;
        for(unsigned int i = 0; i < pointsN3.size(); i++){
            DN3 += coeffsN3[i + 1] * std::exp(x + pointsN3[i] * iter);
        }
        fileOut << std::abs(DN3 / iter - exp0)  <<  ",";
    }
    fileOut << std::endl;

    auto coeffsN4 = calcDerivativeCoef(pointsN4);

    for(auto&& iter : h) {
        double DN4 = coeffsN4[0] * exp0;
        for(unsigned int i = 0; i < pointsN4.size(); i++){
            DN4 += coeffsN4[i + 1] * std::exp(x + pointsN4[i] * iter);
        }
        fileOut << std::abs(DN4 / iter - exp0)  <<  ",";
    }
    fileOut << std::endl;

    auto coeffsN5 = calcDerivativeCoef(pointsN5);

    for(auto&& iter : h) {
        double DN5 = coeffsN5[0] * exp0;
        for(unsigned int i = 0; i < pointsN5.size(); i++){
            DN5 += coeffsN5[i + 1] * std::exp(x + pointsN5[i] * iter);
        }
        fileOut << std::abs(DN5 / iter - exp0)  <<  ",";
    }
    fileOut << std::endl;
    fileOut.close();
}
