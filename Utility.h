#pragma once
#include "Libs.h"
#include "Matrix.h"
#include "CustomException.h"

double getOneVarFuncMinima(std::function<double(double)> f, double a, double b);
Matrix getFunctionGradientInPoint(std::function<double(Matrix)> function, Matrix point);
std::function<double(double)> getDirectionalFunction(std::function<double(Matrix)> function, Matrix dir, Matrix startPoint);
double getBeta(Matrix r_k, Matrix r_k_1);
Matrix getMultiVarFuncLocalMinimum(std::function<double(Matrix)> function, Matrix initPoint);
void getAndShowAllLocalMinimasMultiVarFunc(std::function<double(Matrix)> function, int argsNum, int steps);