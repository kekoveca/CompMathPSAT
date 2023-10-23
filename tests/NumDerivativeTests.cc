#include <gtest/gtest.h>
#include "../src/NumDerivative.hpp"
#include <iostream>
#include <fstream>

const std::string dataPath = "C:/Users/kekoveca/Documents/CompMathPSAT/post_proccessing";

TEST(NumDerivativeTests, main){
    std::ofstream fileOut;
    fileOut.open(dataPath + "/NumDerivativeFolder/errors.txt");

    const unsigned int N = 2;
    const std::array<float, N> points = {1, 2};

    auto coeffs = calcDerivativeCoef(points);

    constexpr double x = 1.;

    constexpr unsigned int N3 = 3;
    constexpr unsigned int N4 = 4;
    constexpr unsigned int N5 = 5;
    
    std::array<double, N3> pointsN3 = {1., 2., 3.};
    std::array<double, N4> pointsN4 = {1., 2., 3., 4.};
    std::array<double, N5> pointsN5 = {1., 2., 3., 4., 5.};

    std::array<double, 16> h{};

    h[0] = 1;

    for(unsigned int i = 0; i < h.size() - 1; i++) {
        h[i + 1] = h[i] * 0.5;
    }

    auto coeffsN3 = calcDerivativeCoef(pointsN3);

    const double exp0 = std::exp(x); 

    for(auto&& iter : h) {
        double DN3 = coeffsN3[0] * exp0;                        // A_0 * f(x_0)
        for(unsigned int i = 0; i < pointsN3.size(); i++){
            DN3 += coeffsN3[i + 1] * pointsN3[i];               // A_i * f(x_0 + ih)
        }
        fileOut<<std::abs(DN3 / iter - exp0)  <<  " ";
    }
    fileOut<<std::endl;

    // DerivativeCoef<double, N4> answN4 = calcDerivativeCoef(pointsN4);
    // //double DN4 = answN4.centralCoef * std::exp(x);

    // for(auto&& it : h){
    //     double DN4 = answN4.centralCoef * std::exp(x);
    //     for(std::uint64_t j = 0; j < pointsN4.size(); j++){
    //         DN4+=answN4.otherCoefs[j] * std::exp(x + pointsN4[j]*it);
    //     }
    //     fileOut<<std::abs(double(DN4)/it - double(std::exp(x)))<<" ";
    // }
    // fileOut<<std::endl;

    // DerivativeCoef<double, N5> answN5 = calcDerivativeCoef(pointsN5);
    // //double DN5 = answN5.centralCoef * std::exp(x);
    // for(auto&& it : h){
    //     double DN5 = answN5.centralCoef * std::exp(x);
    //     for(std::uint64_t j = 0; j < pointsN5.size(); j++){
    //         DN5+=answN5.otherCoefs[j] * std::exp(x + pointsN5[j]*it);
    //     }
    //     fileOut<<std::abs(DN5/it - std::exp(x))<<" ";
    // }
    fileOut<<std::endl;
    fileOut.close();
}
