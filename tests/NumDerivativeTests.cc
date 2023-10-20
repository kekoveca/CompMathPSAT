#include <gtest/gtest.h>
#include "../src/NumDerivative.hpp"


TEST(NumDerivativeTests, FirstTest){
    const unsigned int N = 2;
    const std::array<double, 2> points = {-2, 3};
    DerivativeCoef<double, 2> ans = calcDerivativeCoef(points);
    ASSERT_NO_FATAL_FAILURE();
}