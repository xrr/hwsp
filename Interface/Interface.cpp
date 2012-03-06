#include <iostream>
#include <fstream>
#include <math.h>

#include "Models.h"
#include "Pricers.h"

//concatenate a 1st vector of dates with the vector of swap dates (start & payments)
/*Tree(RateCurve _ratecurve, HullWhite _hullwhite, Swap swap) :  PricerOptions(_ratecurve, _hullwhite) {
Dates blabla;
Tree(_ratecurve, _hullwhite, blabla);
};*/

//std::ofstream(std::ofstream("arbre.csv", std::ios::trunc)) << arbre;
//std::cout << "\n" << arbre;


int main() {

	RateCurve courbe;
	for (int i = 0 ; i<3 ; i++) {
		courbe.times.push_back(2+i);
		courbe.rates.push_back(0.05);
	}
	HullWhite hullwhite(0.05, 0.01);

	PricerGeneric pricer_swap(courbe);
	ClosedFormula pricer_option1(courbe, hullwhite);

	Swap swap(0.05,1,courbe.times);
	Swaption swaption1(swap, true);
	Swaption swaption2(swap, false);


	std::cout << "=== Formule fermee ==="  << "\n";
	std::cout << "Payeur \t\t: " << percentage(pricer_swap.Evaluate(swap)) << "\n";
	std::cout << "Payeuse \t: " << percentage(pricer_option1.Evaluate(swaption1)) << "\n";
	std::cout << "Receveuse \t: " << percentage(pricer_option1.Evaluate(swaption2)) << "\n";

	std::getchar();
	return 0;
}