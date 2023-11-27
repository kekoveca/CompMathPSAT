#include <vector>
#include <type_traits>
#include <algorithm>
#include <map>

const std::vector<double> LEJENDRE_ZEROS_2    = {0.577350269189626};
const std::vector<double> LEJENDRE_WEIGHTS_2  = {1.};

const std::vector<double> LEJENDRE_ZEROS_3    = {0., 0.774596669241483};
const std::vector<double> LEJENDRE_WEIGHTS_3  = {0.888888888888889, 0.55555555555555};

const std::vector<double> LEJENDRE_ZEROS_4    = {0.339981043584856, 0.86113631159405};
const std::vector<double> LEJENDRE_WEIGHTS_4  = {0.652145154862546, 0.34785484513745};

const std::vector<double> LEJENDRE_ZEROS_5    = {0., 0.538469310105683, 0.90617984593866};
const std::vector<double> LEJENDRE_WEIGHTS_5  = {0.568888888888889, 0.478628670499366, 0.236926885056189};

const std::vector<double> LEJENDRE_ZEROS_6    = {0.238619186083197, 0.661209386466265, 0.932469514203152};
const std::vector<double> LEJENDRE_WEIGHTS_6  = {0.467913934572691, 0.360761573048139, 0.171324492379170};

const std::vector<double> LEJENDRE_ZEROS_15   = {0., 0.201194093997435, 0.394151347077563, 0.570972172608539, 
                                                       0.724417731360170, 0.848206583410427, 0.937273392400706, 0.987992518020485};
const std::vector<double> LEJENDRE_WEIGHTS_15 = {0.202578241925561, 0.198431485327111, 0.186161000015562, 0.166269205816994,
                                                       0.139570677926154, 0.107159220467172, 0.070366047488108, 0.030753241996117};

const auto Gauss2 = std::make_pair(LEJENDRE_ZEROS_2, LEJENDRE_WEIGHTS_2);
const auto Gauss3 = std::make_pair(LEJENDRE_ZEROS_3, LEJENDRE_WEIGHTS_3);
const auto Gauss4 = std::make_pair(LEJENDRE_ZEROS_4, LEJENDRE_WEIGHTS_4);
const auto Gauss5 = std::make_pair(LEJENDRE_ZEROS_5, LEJENDRE_WEIGHTS_5);
const auto Gauss6 = std::make_pair(LEJENDRE_ZEROS_6, LEJENDRE_WEIGHTS_6);
const auto Gauss15 = std::make_pair(LEJENDRE_ZEROS_15, LEJENDRE_WEIGHTS_15);


const std::map<std::size_t, std::pair<std::vector<double>, std::vector<double>>> availableQuadratures {
    {2, Gauss2},
    {3, Gauss3},
    {4, Gauss4},
    {5, Gauss5},
    {6, Gauss6},
    {15, Gauss15}
};


template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename T>
using Dif = decltype(std::declval<T>() - std::declval<T>());

/* Функция производит интегрирование на одном отрезке */
template<typename Callable, std::size_t N>
decltype(auto) integrate(const Callable& func, const typename ArgumentGetter<Callable>::Argument& start, const typename ArgumentGetter<Callable>::Argument& end){
    typename ArgumentGetter<Callable>::Argument ans = 0;
    const auto quadrature = availableQuadratures.at(N);
    const auto delta = end - start;
    if(N % 2 != 0) {
        ans += (delta) / 2 * func((end + start) / 2 + quadrature.first[0] * (delta) / 2) * quadrature.second[0];
        for(std::size_t i = 1; i < N / 2 + 1; i++) {
            ans += (delta) / 2 * func((end + start) / 2 + quadrature.first[i] * (delta) / 2) * quadrature.second[i];
            ans += (delta) / 2 * func((end + start) / 2 - quadrature.first[i] * (delta) / 2) * quadrature.second[i];
        }
    }
    else {
        for(std::size_t i = 0; i < N / 2; i++) {
            ans += (delta) / 2 * func((end + start) / 2 + quadrature.first[i] * (delta) / 2) * quadrature.second[i];
            ans += (delta) / 2 * func((end + start) / 2 - quadrature.first[i] * (delta) / 2) * quadrature.second[i];
        }    
    }
    return ans;
}


/* Функция производит интегрирование, разбивая отрезок на подотрезки длиной не более dx */
template<typename Callable, std::size_t N>
decltype(auto) integrate(const Callable& func, const typename ArgumentGetter<Callable>::Argument& start, const typename ArgumentGetter<Callable>::Argument& end, const Dif<typename ArgumentGetter<Callable>::Argument>& dx) {
    typename ArgumentGetter<Callable>::Argument ans = 0;
    const auto quadrature = availableQuadratures.at(N);
    const auto Delta = end - start;
    const unsigned int numberOfIntervals = std::ceil(Delta / dx);
    const Dif<typename ArgumentGetter<Callable>::Argument> exactDx = Delta / numberOfIntervals;
    Dif<typename ArgumentGetter<Callable>::Argument> currStart = start;
    Dif<typename ArgumentGetter<Callable>::Argument> currEnd = start + exactDx;
    for(unsigned int j = 0; j < numberOfIntervals; j++) {
        if(N % 2 != 0) {
            auto delta = currEnd - currStart;
            ans += (delta) / 2 * func((currEnd + currStart) / 2 + quadrature.first[0] * (delta) / 2) * quadrature.second[0];
            for(std::size_t i = 1; i < N / 2 + 1; i++) {
                ans += (delta) / 2 * func((currEnd + currStart) / 2 + quadrature.first[i] * (delta) / 2) * quadrature.second[i];
                ans += (delta) / 2 * func((currEnd + currStart) / 2 - quadrature.first[i] * (delta) / 2) * quadrature.second[i];
            }
            currStart += exactDx;
            currEnd += exactDx;
        }
        else {
            auto delta = currEnd - currStart;
            for(std::size_t i = 0; i < N / 2; i++) {
                ans += (delta) / 2 * func((currEnd + currStart) / 2 + quadrature.first[i] * (delta) / 2) * quadrature.second[i];
                ans += (delta) / 2 * func((currEnd + currStart) / 2 - quadrature.first[i] * (delta) / 2) * quadrature.second[i];
            }
            currStart += exactDx;
            currEnd += exactDx;  
        }
    }
    return ans;
}


template<typename Callable, std::size_t N>
decltype(auto) integrateRichardsonExtrapolate(const Callable& func, 
                                            const typename ArgumentGetter<Callable>::Argument& start, 
                                            const typename ArgumentGetter<Callable>::Argument& end, 
                                            const Dif<typename ArgumentGetter<Callable>::Argument>& dx) {
    decltype(auto) Ih = integrate<double(double), N>(func, start, end, dx);
    decltype(auto) Ih2 = integrate<double(double), N>(func, start, end, dx / 2.);
    return Ih2 + (Ih2 - Ih) / (std::pow(2, 2 * N) - 1);
}


template<typename Callable, std::size_t N>
decltype(auto) integrateRungeRule(const Callable& func, 
                                const typename ArgumentGetter<Callable>::Argument& start, 
                                const typename ArgumentGetter<Callable>::Argument& end, 
                                Dif<typename ArgumentGetter<Callable>::Argument> dx,
                                const typename ArgumentGetter<Callable>::Argument& tolerance) {
    
    decltype(auto) Ih = integrate<double(double), N>(func, start, end, dx);
    decltype(auto) Ih2 = integrate<double(double), N>(func, start, end, dx / 2.);
    const auto theta = 1 / (std::pow(2, 2 * N) - 1);
    decltype(auto) delta = theta * std::abs(Ih - Ih2);
    decltype(Ih) ans;
    while(delta > tolerance) {
        dx = dx / 2.;
        Ih = integrate<double(double), N>(func, start, end, dx);
        Ih2 = integrate<double(double), N>(func, start, end, dx / 2.);
        delta = theta * std::abs(Ih2 - Ih);  
    }
    return Ih2 + delta;
}