#include <gtest/gtest.h>
#include "../src/RungeKutta.hpp"
#include <iostream>
#include <fstream>
#include <array>
#include <Eigen/Dense>

const std::string dataPath = "/Users/kekoveca/Documents/MIPT/CompMath/post_proccessing/RungeKuttaFolder";

TEST(RungeKuttaTests, oscillatorTest) {
    typename Oscillator::Argument initArg = 0;
    typename Oscillator::State initState = {1., 0.};
    typename Oscillator::StateAndArg initStateAndArg = {initState, initArg};
    auto ans = integrate<RK4Table, Oscillator>()
}