#include <gtest/gtest.h>
#include "../src/Spline.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>

const std::string dataPath = "/Users/kekoveca/Documents/MIPT/CompMath/post_proccessing/SplineFolder";

TEST(ThreeDiagonalMatrixTests, mainTest) {
    std::vector<double> mD = {42, 34, 15};
    std::vector<double> lD = {17, 14};
    std::vector<double> uD = {19, 30};
    const ThreeDiagonalMatrix<double> matrix(mD, lD, uD);
    const std::vector<double> d = {18, 73, 43};
    auto ans = solve(matrix, d);
    ASSERT_NEAR(ans[0], -5., 1e-13);
    ASSERT_NEAR(ans[1], 12., 1e-13);
    ASSERT_NEAR(ans[2], -8.33333333333333, 1e-13);
}

TEST(CubicSplineTests, basicInitTest) {
    std::vector<double> points = {0, 1, 2, 4};
    std::vector<double> values = {1, 3, 9, 81};
    CubicSpline<double, double> spline(points, values, 0, 0);
    auto h = spline.getSteps();
    auto splitDiffs1 = spline.getSplitDiffs1();
    auto splitDiffs2 = spline.getSplitDiffs2MultBySix();
    ASSERT_EQ(h[0], 1);
    ASSERT_EQ(h[1], 1);
    ASSERT_EQ(splitDiffs1[0], 2);
    ASSERT_EQ(splitDiffs1[1], 6);
    ASSERT_EQ(splitDiffs2[0], 2 * 6);
}

TEST(CubicSplineTests, mainNaturalTest) {
    std::ofstream fileOutNatural;

    fileOutNatural.open(dataPath + "/errNatural.txt");

    const double leftBound = 0.; 
    const double rightBound = 10.; 
    const unsigned int numberOfExperiments = 6;
    const std::array<unsigned int, numberOfExperiments> numberOfPoints = {5, 10, 20, 40, 80, 160};

    std::array<std::vector<double>, numberOfExperiments> matrixOfPoints;
    std::array<std::vector<double>, numberOfExperiments> matrixOfValues;

    for(unsigned int i = 0; i < numberOfExperiments; i++) {
        matrixOfPoints[i].resize(numberOfPoints[i]);
        matrixOfValues[i].resize(numberOfPoints[i]);
        double multiplier = (rightBound - leftBound) / (numberOfPoints[i] - 1.);
        for(unsigned int j = 0; j < numberOfPoints[i]; j++) {
            matrixOfPoints[i][j] = leftBound + j * multiplier;
            matrixOfValues[i][j] = std::exp(matrixOfPoints[i][j]);
        }
    }

    std::array<double, 1000> pointsForErr;
    std::array<double, 1000> valuesTrue;
    double multiplier = (rightBound - leftBound) / (1000 - 1); 

    for(unsigned int i = 0; i < 1000; i++) {
        pointsForErr[i] = leftBound + multiplier * i;
        valuesTrue[i] = std::exp(pointsForErr[i]);
    }

    for(unsigned int i = 0; i < numberOfExperiments; i++) {
        double error = 0;
        CubicSpline spline(matrixOfPoints[i], matrixOfValues[i], 0., 0.);
        for(unsigned int j = 0; j < 1000; j++) {
            auto interpolatedValue = spline.interpolate(pointsForErr[j]);
            error = std::max(std::abs(interpolatedValue - valuesTrue[j]), error);
        }
        fileOutNatural << error << ",";
    }
    
    fileOutNatural << "\n";
}

TEST(CubicSplineTests, mainGeneralTest) {
    std::ofstream fileOutGeneral;

    fileOutGeneral.open(dataPath + "/errGeneral.txt");

    const double leftBound = 0.; 
    const double rightBound = 10.; 
    const unsigned int numberOfExperiments = 6;
    const std::array<unsigned int, numberOfExperiments> numberOfPoints = {5, 10, 20, 40, 80, 160};

    std::array<std::vector<double>, numberOfExperiments> matrixOfPoints;
    std::array<std::vector<double>, numberOfExperiments> matrixOfValues;

    for(unsigned int i = 0; i < numberOfExperiments; i++) {
        matrixOfPoints[i].resize(numberOfPoints[i]);
        matrixOfValues[i].resize(numberOfPoints[i]);
        double multiplier = (rightBound - leftBound) / (numberOfPoints[i] - 1.);
        for(unsigned int j = 0; j < numberOfPoints[i]; j++) {
            matrixOfPoints[i][j] = leftBound + j * multiplier;
            matrixOfValues[i][j] = std::exp(matrixOfPoints[i][j]);
        }
    }

    std::array<double, 1000> pointsForErr;
    std::array<double, 1000> valuesTrue;
    double multiplier = (rightBound - leftBound) / (1000 - 1); 

    for(unsigned int i = 0; i < 1000; i++) {
        pointsForErr[i] = leftBound + multiplier * i;
        valuesTrue[i] = std::exp(pointsForErr[i]);
    }

    double leftSecondDeriv = 1.;
    double rightSecondDeriv = std::exp(10.);

    for(unsigned int i = 0; i < numberOfExperiments; i++) {
        double error = 0;
        CubicSpline spline(matrixOfPoints[i], matrixOfValues[i], leftSecondDeriv, rightSecondDeriv);
        for(unsigned int j = 0; j < 1000; j++) {
            auto interpolatedValue = spline.interpolate(pointsForErr[j]);
            error = std::max(std::abs(interpolatedValue - valuesTrue[j]), error);
        }
        fileOutGeneral << error << ",";
    }
    
    fileOutGeneral << "\n";
}