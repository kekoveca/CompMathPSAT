#include <gtest/gtest.h>
#include "../src/NonLinearEquations.hpp"
#include <iostream>
#include <fstream>
#include <array>

const std::string dataPath = "/Users/kekoveca/Documents/MIPT/CompMath/post_proccessing/NonLinearEquationsFolder";

TEST(NonLinearEquationsTests, KeplerTest) {
    std::ofstream fileOut;
    fileOut.open(dataPath + "/errKepler.txt");

    std::array<double, 3> eccs = {0.5, 0.2, 0.1};
    auto answer_true = keplerSolver(0.8, M_PI_4, 10, 1e-15);
    unsigned int i_max;
    for(unsigned int i = 1; ; i++) {
        auto answer = keplerSolver(0.8, M_PI_4, i, 1e-15);
        i_max = i;
        if(answer.first) break;
        fileOut << std::abs(answer.second - answer_true.second) << ",";
    }
    fileOut << "\n";
    for(auto&& ecc: eccs) {
        auto answer_true = keplerSolver(ecc, M_PI_4, 10, 1e-15);
        for(unsigned int i = 1; i < i_max; i++) {
            auto answer = keplerSolver(ecc, M_PI_4, i, 1e-15);
            fileOut << std::abs(answer.second - answer_true.second) << ",";
        }
        fileOut << "\n";
    }

    fileOut.close();
}

// подставили y из второго уравнения в первое => F(x) = x^2 + tg^2(x) - 1
double F(double x) {
    double tan = std::tan(x);
    return x * x + tan * tan - 1;
}

TEST(NonLinearEquationsTests, mpiTest) {
    // на отрезке [0, 1] имеем m = 0, M ~ 12 => tau_opt ~ -1/6. Минус из-за того, что f' > 0
    // решения системы - пересечения окружности с центром в начале координат и тангенса, отсюда сразу понятно,
    // что численно нужно искать только один корень

    std::ofstream fileOut;
    fileOut.open(dataPath + "/solutions.txt");

    auto ans = solve<double(double), double>(F, -1./6, 0.8, 100, 1e-6);
    auto x1 = ans.second;
    auto y1 = std::tan(x1);
    auto x2 = - x1; 
    auto y2 = - y1; 
    fileOut << "Решения:\n(" << x1 << ", " << y1 <<")\n(" << x2 << ", " << y2 << ")\n";

    fileOut.close();
}