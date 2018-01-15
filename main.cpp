#include <iostream>
#include <cstdio>
#include "Libs.h"
#include "Utility.h"


int main(int argc, char* argv[]) {
	std::function<double(Matrix)> func = [](Matrix point) {
		return pow(2 - point(0, 0), 2) + 47.*pow(point(1, 0) - pow(point(0, 0), 2), 2);
	};

	Matrix initPoint(2, 1);
	initPoint(0, 0) = 1;
	initPoint(1, 0) = 1;

	try {
		getAndShowAllLocalMinimasMultiVarFunc(func, 2, 20);
	}
	catch (std::exception e) {
		std::cout << e.what() << std::endl;
	}
	

	_getch();
	return 0;
}