#include <iostream>
#include <fstream>
#include <math.h>

#include "Models.h"
#include "Pricers.h"

int main() {

	//for (int i = 0 ; i<10 ; i++) {
	//	courbe.times.push_back(0.02*(1+i));
	//	//courbe.rates.push_back(0.08-0.05*exp(-0.18*(0.2*(2+i))));
	//courbe.rates.push_back(0.05);
	//}

	/*std::cout << courbe.zerocoupon(0) << std::endl;
	std::cout << courbe.zerocoupon(0.5) << std::endl;
	std::cout << courbe.zerocoupon(1) << std::endl;
	std::cout << courbe.zerocoupon(1.5) << std::endl;
	std::cout << courbe.zerocoupon(2) << std::endl;
	std::cout << courbe.zerocoupon(3) << std::endl;
	double a = 8;*/

	//========= Tests fonction rate(t) ======
	/*for (int i = 0 ; i<2 ; i++) {
		courbe.times.push_back(1+i);
		courbe.rates.push_back(2*(1+i));
	}*/
	//std::cout << courbe.rate(0.5) << std::endl; -> 2 expected
	//std::cout << courbe.rate(1) << std::endl; -> 2 expected
	//std::cout << courbe.rate(1.5) << std::endl; -> 2 expected
	//std::cout << courbe.rate(2) << std::endl; -> 4 expected
	//std::cout << courbe.rate(3) << std::endl; -> 4 expected
	//end test fonction rate(t)
	//========================================


	RateCurve courbe;

	for (int i = 0 ; i<10 ; i++) {
		courbe.times.push_back(0.02*(1+i));
		//courbe.rates.push_back(0.08-0.05*exp(-0.18*(0.2*(2+i))));
	courbe.rates.push_back(0.05);
	}

	Dates d(courbe.times);
	Swap swap(0.05,0.01,d);
	Swaption swaption1(swap, d, true);
	Swaption swaption2(swap, d, false);
	HullWhite hullwhite(0.1, 0.01);
	
	std::cout << "Prix du Swaption payer :" << ClosedFormula (courbe, hullwhite).Evaluate(swaption1) << "\n";
	std::cout << "Prix du Swaption receiver :" << ClosedFormula (courbe, hullwhite).Evaluate(swaption2) << "\n";

	std::cout << "Prix du Swap Formule fermee :" << PricerGeneric (courbe).Evaluate(swap) << "\n";
	//std::cout << "Prix Zero Coupon en 0:" << courbe.zerocoupon(0) << "\n";

	////std::ofstream(std::ofstream("arbre.csv", std::ios::trunc)) << arbre;
	////std::cout << "\n" << arbre;

	std::getchar();
	return 0;
}