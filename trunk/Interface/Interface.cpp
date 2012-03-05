#include <iostream>
#include <fstream>
#include <math.h>

#include "Models.h"
#include "Pricers.h"

int main() {

	RateCurve courbe;

	for (int i = 0 ; i<2 ; i++) {
		courbe.times.push_back(1*(2+i));
		courbe.rates.push_back(0.08-0.05*exp(-0.18*(1*(2+i))));
	}

	Dates d(courbe.times);
	Swap swap(0.07,1,d);
	//RateCurve courb2;
	//	courbe2.times.push_back(0);
	//Dates mat_swaption(
	Swaption swaption1(swap, d, true);
	Swaption swaption2(swap, d, false);
	HullWhite hullwhite(0.1, 0.01);
	
	std::cout << "Prix du Swaption payer :" << ClosedFormula (courbe, hullwhite).Evaluate(swaption1) << "\n";
	std::cout << "Prix du Swaption receiver :" << ClosedFormula (courbe, hullwhite).Evaluate(swaption2) << "\n";

	std::cout << "Prix du Swap Formule fermee :" << PricerGeneric (courbe).Evaluate(swap) << "\n";
	
	//std::ofstream(std::ofstream("arbre.csv", std::ios::trunc)) << arbre;
	//std::cout << "\n" << arbre;

	std::getchar();
	return 0;
}