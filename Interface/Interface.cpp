#include <iostream>
#include <fstream>
#include <math.h>

#include "Models.h"
#include "Pricers.h"

int main() {

	RateCurve courbe;

	for (int i = 0 ; i<10 ; i++) {
		courbe.times.push_back(0.2*(1+i));
		courbe.rates.push_back(0.08-0.05*exp(-0.18*(0.2*(1+i))));
	}

	Dates d(courbe.times);
	Swap swap(0.05, d);
	Swaption swaption(swap, d, true);
	HullWhite hullwhite(0.1, 0.01);
	
	std::cout << ClosedFormula (courbe, hullwhite).Evaluate(swaption);
		
	
	//std::ofstream(std::ofstream("arbre.csv", std::ios::trunc)) << arbre;
	//std::cout << "\n" << arbre;

	std::getchar();
	return 0;
}