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

	for (int i = 0 ; i<3 ; i++) {
		//courbe.times.push_back(0.02*(1+i));
		//courbe.rates.push_back(0.08-0.05*exp(-0.18*(0.2*(2+i))));
		//courbe.rates.push_back(0.05);
		courbe.times.push_back(2+i);
		courbe.rates.push_back(0.05);
	
	}
	Swap swap(0.05,1,courbe.times);
	Swaption swaption1(swap, true);
	Swaption swaption2(swap, false);
	HullWhite hullwhite(0.05, 0.01);
	
	//Dates d(courbe.times);
	//Swap swap(0.05,0.4,d);

	
	//concatenate a 1st vector of dates with the vector of swap dates (start & payments)
	/*Tree(RateCurve _ratecurve, HullWhite _hullwhite, Swap swap) :  PricerOptions(_ratecurve, _hullwhite) {
		Dates blabla;
		Tree(_ratecurve, _hullwhite, blabla);
	};*/

	std::cout << "Prix du Swap Formule fermee :" << PricerGeneric (courbe).Evaluate(swap) << "\n";
	std::cout << "Prix du Swaption Formule fermee payer :" << ClosedFormula (courbe, hullwhite).Evaluate(swaption1) << "\n";
	std::cout << "Prix du Swaption Formule fermee receiver :" << ClosedFormula (courbe, hullwhite).Evaluate(swaption2) << "\n";
	//std::cout << "Prix Zero Coupon en 0:" << courbe.zerocoupon(0) << "\n";

	////std::ofstream(std::ofstream("arbre.csv", std::ios::trunc)) << arbre;
	////std::cout << "\n" << arbre;

	std::getchar();
	return 0;
}