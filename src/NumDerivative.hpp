#include <iostream>
#include <Eigen/Dense>
#include <cmath>

template<typename RealType, std::size_t N>
struct DerivativeCoef {
    RealType centralCoef;
    std::array<RealType, N> otherCoefs;
};

template<typename RealType, std::size_t N>
// DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points) noexcept {
Eigen::VectorXf calcDerivativeCoef(const std::array<RealType, N>& points) noexcept {  
    Eigen::Vector<RealType, N + 1> b;
    b.setZero();
    b(1) = 1;

    Eigen::Matrix<RealType, N + 1, N + 1> A;
    Eigen::Array<RealType, 1, N + 1> tmp;
    tmp.setOnes();
    A.row(0) = tmp;

    auto pointsCopy = points;

    for(unsigned int j = 1; j < N + 1; j++) {
        A(j, 0) = 0;
        for(unsigned int i = 1; i < N + 1; i++) {
            A(j, i) = std::pow(pointsCopy[i - 1], j);
        }
    }

    typedef Matrix<RealType, N, Dynamic

    Eigen::ColPivHouseholderQR<Eigen::MatrixXf> dec(A);
    Eigen::VectorXf allCoefs = dec.solve(b);

    DerivativeCoef<RealType, N> ans;

    for(unsigned int i = 0; i < 1; i++) {
        ans.otherCoefs[i] = allCoefs[i];
    }

    ans.centralCoef = allCoefs(1);
    
    for(unsigned int i = 2; i < N + 1; i++) {
        ans.otherCoefs[i - 1] = allCoefs[i];
    }

    return {};
}