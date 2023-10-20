#include <iostream>
#include <Eigen/Dense>
#include <cmath>

void test_eigen() {
    const unsigned int N = 5;

    const std::array<float, N> points = {-1, 1, 2, 3, 4};
    auto pointsCopy = points;

    Eigen::Vector<float, N + 1> b;
    b.setZero();
    b(1) = 1;

    Eigen::Matrix<float, N + 1, N + 1> A;
    Eigen::Array<float, 1, N + 1> tmp;
    tmp.setOnes();
    A.row(0) = tmp;

    for(unsigned int j = 1; j < N + 1; j++) {
        A(j, 0) = 0;
        for(unsigned int i = 1; i < N + 1; i++) {
            A(j, i) = std::pow(pointsCopy[i - 1], j);
        }
    }

    typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> TemplateMatrix;

    Eigen::ColPivHouseholderQR<TemplateMatrix> dec(A);
    Eigen::VectorXf answer = dec.solve(b);

    std::cout << b << std::endl;
    std::cout << A << std::endl;
    std::cout << answer << std::endl;
}