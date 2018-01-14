#include <iostream>
#include <cstdio>
#include "Libs.h"
#include "Utility.h"


/*
IDEA
Funkcja f(x) na przedziale [a,b] posiada minimum.
Poprzez kolejne podziały przedziału [a,b] otrzymamy te minimum.
Szukamy wartorsci xL i xR, takich, że: a<xL<xR<b
Jeśli f(xL) > f(xR), wówczas minimum znajduje się na przedziale [xL, b]
jeśli nie, to minimum znajduje się na przedziale [a, xR]
*/
double getOneVarFuncMinima(std::function<double(double)> f, double a, double b) {

	//Wspołczynnik złotego podziału
	double gsFactor = (sqrt(5) + 1.) / 2.;

	//Tolerancja
	double epsilon = 1e-6;

	//c - lewa próbka
	//d - prawa próbka
	double xL = b - (b - a) / gsFactor;
	double xR = a + (b - a) / gsFactor;

	//Ponieważ minimum będzie znajdowało się pomiędzy xL i xR, to chcemy, aby wartosc bezwzgledna z ich roznicy byla mniejsza niz epsilon
	while (std::abs(xL - xR) > epsilon) {
		if (f(xL) < f(xR))
			b = xR; //wybieramy przedział [a,xR]
		else
			a = xL; //wybieramy przedział [xL,b]

		xL = b - (b - a) / gsFactor;
		xR = a + (b - a) / gsFactor;
	}
	return (a + b) / 2; //zwaracamy wartość środkową dla zaakceptowanego przedziału
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
	std::function<double(double)> tFunction = [=](double alfa) mutable  {
		return function(startPoint + dir.getMultipliedByScalar(alfa));
	};
	return tFunction;
}

double getBeta(Matrix r_k, Matrix r_k_1) {
	return (r_k_1.getScalarProduct(r_k_1 - r_k)) / (r_k.getScalarProduct(r_k));
	//return (r_k_1.getScalarProduct(r_k_1)) / (r_k.getScalarProduct(r_k));
}

Matrix getLocalMinimum(std::function<double(Matrix)> function, Matrix initPoint) {
	
	//Preparation
	Matrix x_k(initPoint); //copy constructor
	Matrix r_k(initPoint.getSize());
	Matrix p_k(initPoint.getSize());


	r_k = getFunctionGradientInPoint(function, x_k).getMultipliedByScalar(-1); //r_i set to -gradient
	p_k = r_k; 
	

	//Step 1:
	auto dirFunc = getDirectionalFunction(function, p_k, initPoint); //directional function
	double dirMinima = getOneVarFuncMinima(dirFunc, -200, 200); 
	Matrix x_k_1 = x_k + p_k.getMultipliedByScalar(dirMinima);

	
	//Step 2:
	Matrix r_k_1 = getFunctionGradientInPoint(function, x_k_1).getMultipliedByScalar(-1);
	
	//Step 3:
	double beta = getBeta(r_k, r_k_1);
	
	//Step 4: getting next direction
	Matrix p_k_1 = r_k_1 + p_k.getMultipliedByScalar(beta);
	

	
	
	
	for (int i = 0; i < 40; i++) {
		x_k = x_k_1;
		r_k = r_k_1;
		p_k = p_k_1;
		dirFunc = getDirectionalFunction(function, p_k, x_k); //directional function
		dirMinima = getOneVarFuncMinima(dirFunc, -1000, 1000); 
		x_k_1 = x_k + p_k.getMultipliedByScalar(dirMinima);
		r_k_1 = getFunctionGradientInPoint(function, x_k_1).getMultipliedByScalar(-1);
		beta = getBeta(r_k, r_k_1);
		p_k_1 = r_k_1 + p_k.getMultipliedByScalar(beta);
	}
	return x_k;
}


int main(int argc, char* argv[]) {

	//F(x,y) = (2-x)^2 +47(y-x^2)^2
	//(x + 2)2 + (y + 3)2
	std::function<double(Matrix)> func = [](Matrix point) {
		//return pow((point(0, 0) + 2), 2) + pow((point(1, 0) + 3), 2);
		return pow(2 - point(0, 0), 2) + 47.*pow(point(1, 0) - pow(point(0, 0), 2), 2);
	};

	Matrix initPoint(2, 1);
	initPoint(0, 0) = 1;
	initPoint(1, 0) = 1;

	
	try {
		std::cout << "Minimum = " << getLocalMinimum(func, initPoint);
	}
	catch (std::exception e) {
		std::cout << e.what() << std::endl;
	}
	

	_getch();
	return 0;
}