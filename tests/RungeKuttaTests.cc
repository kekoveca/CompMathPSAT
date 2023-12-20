#include <gtest/gtest.h>
#include "../src/RungeKutta.hpp"
#include <iostream>
#include <fstream>
#include <array>

const std::string dataPath = "C:\\Users\\kekoveca\\Documents\\CompMathPSAT\\post_proccessing\\RungeKuttaFolder";


TEST(RungeKuttaTests, mainFirstTest){
    std::ofstream errFile;
    errFile.open(dataPath + "/err.txt");
    double endTime = 5.;

    auto func = [](auto x){return std::pow(x, 4) / 4.;};
    double step = 2;
    for(std::size_t i = 0; i <= 12; i++){
        step = step * 0.5;
        int stepNumber = endTime / step;
        std::vector<double> resultTrue;
        double err = 0;
        for(std::size_t j = 0; j <= stepNumber; j++){
            resultTrue.push_back(func(step * j));
        }
        resultTrue.push_back(func(endTime));
        auto result = integrate<RK4Table, tCube>(tCube::StateAndArg{tCube::State {0.}, 0}, endTime, step, tCube{});
        for(std::size_t j = 0; j < result.size(); j++){

            err = std::max(err, std::abs((result[j].state(0) - resultTrue[j])));
        }
        errFile << err << ",";
    }
    errFile.close();
}

TEST(RungeKuttaTests, mainSecondTest){
    std::ofstream errFile;
    errFile.open(dataPath + "/err1.txt");
    double endTime = 5.;

    auto func = [](auto x){return std::sin(x);};
    double step = 2;
    for(std::size_t i = 0; i <= 12; i++){
        step = step * 0.5;
        int stepNumber = endTime / step;
        std::vector<double> resultTrue;
        double err = 0;
        for(std::size_t j = 0; j <= stepNumber; j++){
            resultTrue.push_back(func(step * j));
        }

        resultTrue.push_back(func(endTime));

        auto result = integrate<RK4Table, Oscillator>(Oscillator::StateAndArg{Oscillator::State{0., 1.}, 0.}, endTime, step, Oscillator{});

        for(std::size_t j = 0; j < result.size(); j++){
            err = std::max(err, std::abs((result[j].state(0) - resultTrue[j])));
        }
        errFile << err << ",";
    }
    errFile.close();
}