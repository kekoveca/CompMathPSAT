#include <iostream>
#include <Eigen/Dense>
#include <cmath>

template<typename RealType, unsigned int N>
struct DerivativeCoef {
    RealType centralCoef;
    std::array<RealType, N> otherCoefs;
};

template<typename RealType, unsigned int N>
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points) noexcept {
    const unsigned int numberOfEquations = points.size() + 1;
    Eigen::Array<RealType, numberOfEquations, 1> b;
    b.setZero();
    b(1) = 1;

    Eigen::Matrix<RealType, numberOfEquations, numberOfEquations> A;
    Eigen::Array<RealType, 1, numberOfEquations> tmp;
    tmp.setOnes();
    A(0) = tmp;

    auto pointsCopy = points;

    for(unsigned int j = 1; j < numberOfEquations; j++) {
        A(j, 0) = 0;
        for(unsigned int i = 1; i < numberOfEquations; i++) {
            A(j, i) = std::pow(pointsCopy[i - 1], j);
        }
    }

    std::array<RealType, N + 1> allCoefs = A.colPivHouseholderQr().solve(b);

    DerivativeCoef<RealType, N> ans;
    ans.centralCoef = allCoefs[1];
    for(unsigned int i = 0; i != 1, i < N + 1; i++) {
        ans.otherCoefs[i] = allCoefs[i];
    }

    return ans;
}