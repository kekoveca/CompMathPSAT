#include <array>
#include <Eigen/Dense>
#include <Eigen/Core>

template<typename xType, typename yType, unsigned int N>
class NewtonInterpolator {
    std::array<yType, N> splitDifferencies;
    std::array<xType, N> points;

    public:
    NewtonInterpolator(const std::array<xType, N> &points, const std::array<yType, N>& values) noexcept: splitDifferencies(values), points(points){
        for(unsigned int i = 0; i < N - 1 ; i++){
            for(int j = N - 2 - i; j >= 0; j--){
                yType splitDifferenciesTemp = (splitDifferencies[j + i + 1] - splitDifferencies[j + i])/(points[j + i + 1] - points[j]);
                splitDifferencies[j + i + 1] = splitDifferenciesTemp;
            }
        }
    }

    std::array<xType, N> getPoints() const{
        return points;
    }

    std::array<yType, N> getsplitDifferencies() const{
        return splitDifferencies;
    }

    yType interpolate(const xType& x) const noexcept{
        yType f = splitDifferencies.back();
        for(int i = N - 1; i > 0;  i--){
            f = splitDifferencies[i - 1] + (x - points[i - 1]) * f;
        }
        return f;
    }
};


template<typename xType, typename yType, unsigned int N>
class EigenNewtonInterpolator {
    Eigen::Array<xType, N, 1> splitDifferencies;
    Eigen::Array<yType, N, 1> points;


    public:
    EigenNewtonInterpolator(const Eigen::Array<xType, N, 1> &points, const Eigen::Array<yType, N, 1>& values) noexcept: splitDifferencies(values), points(points){
        for(unsigned int i = 0; i < N - 1 ; i++){
            for(int j = N - 2 - i; j >= 0; j--){
                yType splitDifferenciesTemp = (splitDifferencies[j + i + 1] - splitDifferencies[j + i])/(points[j + i + 1] - points[j]);
                splitDifferencies[j + i + 1] = splitDifferenciesTemp;
            }
        }
    }

    Eigen::Array<xType, N, 1> getPoints() const{
        return points;
    }

    Eigen::Array<yType, N, 1>  getsplitDifferencies() const{
        return splitDifferencies;
    }

    yType interpolate(const xType& x) const noexcept{
        yType f = splitDifferencies.tail(1)(0);
        for(int i = N - 1; i > 0;  i--){
            f = splitDifferencies[i - 1] + (x - points[i - 1]) * f;
        }
        return f;
    }
};