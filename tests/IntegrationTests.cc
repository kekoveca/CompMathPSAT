#include <gtest/gtest.h>
#include "../src/Integration.hpp"
#include <iostream>
#include <fstream>

const std::string dataPath = "/Users/kekoveca/Documents/MIPT/CompMath/post_proccessing/IntegrationFolder";

TEST(IntegrationTests, basicTest) {
    double answer = 1. - std::cos(double(10.));
    constexpr double start = 0;
    constexpr double end = 10;
    constexpr std::array<std::size_t, 6> nums = {2, 3, 4, 5, 6, 15};
    std::cout << std::abs(answer - integrate<double(double), 2>(std::sin, start, end)) << "\n";
    std::cout << std::abs(answer - integrate<double(double), 3>(std::sin, start, end)) << "\n";
    std::cout << std::abs(answer - integrate<double(double), 4>(std::sin, start, end)) << "\n";
    std::cout << std::abs(answer - integrate<double(double), 5>(std::sin, start, end)) << "\n";
    std::cout << std::abs(answer - integrate<double(double), 6>(std::sin, start, end)) << "\n";
    std::cout << std::abs(answer - integrate<double(double), 15>(std::sin, start, end)) << "\n";
}

TEST(IntegrationTests, firstTest){
    std::ofstream fileOut;
    fileOut.open(dataPath + "/err.txt");
    constexpr double start = 0;
    constexpr double end = 10;
    constexpr std::array dx{8., 4., 2., 1., 1./2., 1./4., 1./8., 1./16., 1./32., 1./64., 1./128., 1./256., 1./512., 1./1024.};
    double answer = 1. - std::cos(double(10.));
    for(auto&& it: dx){
        fileOut << std::abs(answer - integrate<double(double), 2>(std::sin, start, end, it)) << ",";
    }
    fileOut << "\n";
    for(auto&& it: dx){
        fileOut << std::abs(answer - integrate<double(double), 3>(std::sin, start, end, it)) << ",";
    }
    fileOut << "\n";
    for(auto&& it: dx){
        fileOut << std::abs(answer - integrate<double(double), 4>(std::sin, start, end, it)) << ",";
    }
    fileOut << "\n";
    for(auto&& it: dx){
        fileOut << std::abs(answer - integrate<double(double), 5>(std::sin, start, end, it)) << ",";
    }
    fileOut << "\n";
    for(auto&& it: dx){
        fileOut << std::abs(answer - integrate<double(double), 6>(std::sin, start, end, it)) << ",";
    }
    fileOut << "\n";
    for(auto&& it: dx){
        fileOut << std::abs(answer - integrate<double(double), 15>(std::sin, start, end, it)) << ",";
    }
    fileOut.close();
}