#include <gtest/gtest.h>
#include "../src/NumDerivative.hpp"
#include "../src/test.hpp"


TEST(NumDerivativeTests, FirstTest){
    const long N = 2;
    const std::array<long, 2> points = {-1, 1};
    // calcDerivativeCoef(points);
    test_eigen();
}