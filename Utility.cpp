#include "Utility.h"

double getOneVarFuncMinima(std::function<double(double)> f, double a, double b) {

	//Golden section factor
	double gsFactor = (sqrt(5) + 1.) / 2.;

	//Tolerance
	double epsilon = 1e-6;

	//xL - left sample
	//xR - right sample
	double xL = b - (b - a) / gsFactor;
	double xR = a + (b - a) / gsFactor;

	//Our minimum will be between xL and xR, so we want to the absolute value of their difference to be less than epsilon
	while (std::abs(xL - xR) > epsilon) {
		if (f(xL) < f(xR))
			b = xR; //choose [a,xR] interval
		else
			a = xL; //choose [xL,b] interval

		xL = b - (b - a) / gsFactor;
		xR = a + (b - a) / gsFactor;
	}
	return (a + b) / 2; //return average value of accepted interval
}
Matrix getFunctionGradientInPoint(std::function<double(Matrix)> function, Matrix point) {

	Matrix gradient(point.getSize());
	double h = sqrt(DBL_EPSILON); //sqrt of smallest value in cpp; it's recommended to take square

	Matrix diffVarFactor(point.getSize());
	for (int i = 0; i < gradient.getSize().getRowsNum(); i++) {
		diffVarFactor.init();
		diffVarFactor(i, 0) = h;

		gradient(i, 0) = (function(point + diffVarFactor) - function(point - diffVarFactor)) / (2.*h); //partial derivative in point
	}
	return gradient;
}
std::function<double(double)> getDirectionalFunction(std::function<double(Matrix)> function, Matrix dir, Matrix startPoint) {
	std::function<double(double)> tFunction = [=](double alfa) mutable {
		return function(startPoint + dir.getMultipliedByScalar(alfa));
	};
	return tFunction;
}
double getBeta(Matrix r_k, Matrix r_k_1) {
	return (r_k_1.getScalarProduct(r_k_1 - r_k)) / (r_k.getScalarProduct(r_k));
}
Matrix getMultiVarFuncLocalMinimum(std::function<double(Matrix)> function, Matrix initPoint) {

	//Preparation
	Matrix x_k(initPoint); //copy constructor
	Matrix r_k(initPoint.getSize());
	Matrix p_k(initPoint.getSize());

	r_k = getFunctionGradientInPoint(function, x_k).getMultipliedByScalar(-1); //r_i set to -gradient
	p_k = r_k;

	//Step 1:
	auto dirFunc = getDirectionalFunction(function, p_k, initPoint); //directional function
	double dirMinima = getOneVarFuncMinima(dirFunc, -1000, 1000);
	Matrix x_k_1 = x_k + p_k.getMultipliedByScalar(dirMinima);

	//Step 2:
	Matrix r_k_1 = getFunctionGradientInPoint(function, x_k_1).getMultipliedByScalar(-1);

	//Step 3:
	double beta = getBeta(r_k, r_k_1);

	//Step 4: getting next direction
	Matrix p_k_1 = r_k_1 + p_k.getMultipliedByScalar(beta);

	//for 200 iterations we will have a good precision of seeked minimum
	for (int i = 0; i < 200; i++) {
		x_k = x_k_1;
		r_k = r_k_1;
		p_k = p_k_1;
		dirFunc = getDirectionalFunction(function, p_k, x_k); //directional function
		dirMinima = getOneVarFuncMinima(dirFunc, -2000, 2000);
		x_k_1 = x_k + p_k.getMultipliedByScalar(dirMinima);
		r_k_1 = getFunctionGradientInPoint(function, x_k_1).getMultipliedByScalar(-1);
		beta = getBeta(r_k, r_k_1);
		p_k_1 = r_k_1 + p_k.getMultipliedByScalar(beta);
	}
	return x_k;
}
void getAndShowAllLocalMinimasMultiVarFunc(std::function<double(Matrix)> function, int argsNum, int steps) {


	std::vector<Matrix> minimas;
	Matrix initVec(argsNum, 1);

	for (int i = 0; i < steps; i++) {
		initVec.initRandomValues(-20, 20);
		Matrix currLocalMinimum = getMultiVarFuncLocalMinimum(function, initVec);
		if (!(std::find(minimas.begin(), minimas.end(), currLocalMinimum) != minimas.end())) {
			bool foundNan = false;
			for (int i = 0; i < currLocalMinimum.getSize().getRowsNum(); i++) {
				if (std::isnan(currLocalMinimum(i, 0))) {
					foundNan = true;
					break;
				}
			}
			if(!foundNan)
				minimas.push_back(currLocalMinimum);
		}
			
	}

	std::cout << "Local minimas of function: " << std::endl;
	for (int i = 0; i < minimas.size(); i++)
		std::cout << i+1 << ") " << std::endl << minimas.at(i) << std::endl;
}
