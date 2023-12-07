#include <Eigen/Dense>
#include <array>
#include <vector>

/* Это таблица Бутчера для метода Рунге-Кутты 4 порядка. Я ее не заполнил */
struct RK4Table{
    static constexpr unsigned int stages = 4;
    static constexpr std::array<std::array<double, stages>, stages> table = {{{0., 0., 0., 0.}, 
                                                                        {0.5, 0., 0., 0.}, 
                                                                        {0., 0.5, 0., 0.}, 
                                                                        {0., 0., 1., 0.}}};
    static constexpr std::array<double, stages> cColumn = {0, 0.5, 0.5, 1.};
    static constexpr std::array<double, stages> bString = {1./6, 1./3, 1./3, 1./6};
};

class Oscillator {

public:

    static constexpr unsigned int dim = 2;  // размерность задачи

    using Argument = double;  // тип аргумента, тип t

    using State = Eigen::Vector<double, dim>;  // состояние

    struct StateAndArg{
        State state;
        Argument arg;
    };

    /*** Вычисляет правую часть ДУ - функцию f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        return Eigen::Vector<double, dim>{stateAndArg.state(1), -stateAndArg.state(0)};
    } 
};

template<typename Table, typename RHS>  // таблица бутчера и класс правой части f
std::vector<typename RHS::StateAndArg> integrate(
    const typename RHS::StateAndArg& initialState,
    const typename RHS::Argument& endTime,
    double step,
    const RHS& rhs
) {
    static constexpr unsigned int dim = rhs.dim;
    rhs.StateAndArg = initialState;
    auto S = Table::stages;
    std::array<Eigen::Vector<double, dim>, S> arrayOfK;
    const Eigen::Vector<double, dim> k1 = rhs.calc(initialState);
    typename RHS::StateAndArg currentTimeState = initialState;

    while(currentTimeState.arg < endTime) {
        arrayOfK[0] = k1;
        auto currentState = currentTimeState;
        for(unsigned int i = 1; i < S; i++) {
            double sumAK = 0;
            for(unsigned int j = 0; j < i - 1; j++) {
                sumAK += Table::table[i][j] * arrayOfK[j];
            }
            currentState = {currentTimeState.state + step * sumAK, currentTimeState.arg + step * Table::cColumn[i]};
            arrayOfK[i] = rhs.calc(currentState);
        }
        double sumBK = 0;
        for(unsigned int i = 0; i < S; i++) {
            sumBK += Table::bString[i] * arrayOfK[i];
        }
        currentTimeState.state = currentTimeState.state + step * sumBK;
        currentTimeState.arg += step;
    }

    return currentTimeState;
}