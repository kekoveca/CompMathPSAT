#include <array>
#include <type_traits>

constexpr std::array<double, 1> LEJENDRE_ZEROS_2    = {0.577350269189626};
constexpr std::array<double, 1> LEJENDRE_WEIGHTS_2  = {1.};

constexpr std::array<double, 2> LEJENDRE_ZEROS_3    = {0., 0.774596669241483};
constexpr std::array<double, 2> LEJENDRE_WEIHGTS_3  = {0.888888888888889, 0.55555555555555};

constexpr std::array<double, 2> LEJENDRE_ZEROS_4    = {0.339981043584856, 0.86113631159405};
constexpr std::array<double, 2> LEJENDRE_WEIHGTS_4  = {0.652145154862546, 0.34785484513745};

constexpr std::array<double, 3> LEJENDRE_ZEROS_5    = {0., 0.538469310105683, 0.90617984593866};
constexpr std::array<double, 3> LEJENDRE_WEIHGTS_5  = {0.568888888888889, 0.478628670499366, 0.236926885056189};

constexpr std::array<double, 3> LEJENDRE_ZEROS_6    = {0.238619186083197, 0.661209386466265, 0.932469514203152};
constexpr std::array<double, 3> LEJENDRE_WEIHGTS_6  = {0.467913934572691, 0.360761573048139, 0.171324492379170};

constexpr std::array<double, 8> LEJENDRE_ZEROS_15   = {0., 0.201194093997435, 0.394151347077563, 0.570972172608539, 
                                                       0.724417731360170, 0.848206583410427, 0.937273392400706, 0.987992518020485};
constexpr std::array<double, 8> LEJENDRE_WEIHGTS_15 = {0.202578241925561, 0.198431485327111, 0.186161000015562, 0.166269205816994,
                                                       0.139570677926154, 0.107159220467172, 0.070366047488108, 0.030753241996117};

const std::array<std::size_t, 6> availableQuadratures= {2, 3, 4, 5, 6, 15};

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

}




// /* Функция производит интегрирование, разбивая отрезок на подотрезки длиной не более dx */
// template<typename Callable, std::size_t N>
// decltype(auto) integrate(   
//     const Callable& func,  // Интегрируемая функция
//     const typename ArgumentGetter<Callable>::Argument& start,  // начало отрезка
//     const typename ArgumentGetter<Callable>::Argument& end,  // конец отрезка
//     const Dif<typename ArgumentGetter<Callable>::Argument>& dx  // Длина подотрезка
//                         );