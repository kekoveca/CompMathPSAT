#include <Eigen/Dense>
#include <cmath>


template<typename RealType, std::size_t N>
auto calcDerivativeCoef(const std::array<RealType, N>& points) {

    Eigen::Matrix<RealType, N + 1, 1> b;
    b.setZero();
    b(1) = 1;

    Eigen::Matrix<RealType, N + 1, N + 1> A;
    Eigen::Matrix<RealType, 1, N + 1> tmp;
    tmp.setOnes();
    A.row(0) = tmp;

    for(unsigned int j = 1; j < N + 1; j++) {
        A(j, 0) = 0;
        for(unsigned int i = 1; i < N + 1; i++) {
            A(j, i) = std::pow(points[i - 1], j);
        }
    }

    typedef Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> RealMatrix;
    typedef Eigen::Vector<RealType, N + 1> AnswArray;

    Eigen::ColPivHouseholderQR<RealMatrix> dec(A);
    AnswArray answerTemporary = dec.solve(b);

    std::array<RealType, N + 1> answer;
    
    for(unsigned int i = 0; i < N + 1; i++) {
        answer[i] = answerTemporary[i];
    }

    return answer;
}