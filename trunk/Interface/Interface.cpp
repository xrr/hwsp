#include <iostream>
#include <fstream>

#include "Node.h"
#include <math.h>

int main() {

	Curve courbe;

	for (int i = 0 ; i<100 ; i++) {
		courbe.times.push_back(0.2*(1+i));
		courbe.rates.push_back(0.08-0.05*exp(-0.18*(0.2*(1+i))));
	}

	Tree arbre(HullWhite(0.1, 0.01, courbe));
	arbre.construct();
	
	std::ofstream(std::ofstream("arbre.csv", std::ios::trunc)) << arbre;
	
	std::cout << "\n" << arbre;

	
	std::getchar();
}

