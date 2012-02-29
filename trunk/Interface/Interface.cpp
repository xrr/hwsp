#include <iostream>
#include "Node.h"

#include <math.h>

int main() {

	std::vector<double> dt;
	for (int i = 0 ; i<3 ; i++) dt.push_back(1);

	HullWhite hw;
	hw.a = 0.1;
	hw.sigma = 0.01;

	Tree arbre(dt,hw);
	arbre.construct();
	std::cout << arbre;
	std::getchar();
}