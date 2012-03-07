#include <iostream>
#include <fstream>
#include <vector>

#include <math.h>

#include "Pricers.h"
#include "Tools.h"


int main () {

	RateCurve courbe;
	for (int i = 0 ; i<3 ; i++) {
		courbe.times.push_back(2+i);
		courbe.rates.push_back(0.05);
	}

	std::vector<double> tree_dates;
	for (int i = 0 ; i<4 ; i++) tree_dates.push_back(0.2*(1+i));
	tree_dates.insert(tree_dates.end(), courbe.times.begin(), courbe.times.end());

	HullWhite hullwhite(0.1, 0.01);


	PricerGeneric pricer_swap(courbe);
	ClosedFormula pricer_option1(courbe, hullwhite);
	Tree pricer_option2(courbe,hullwhite, tree_dates);


	//Swap swap(0.05,1,courbe.times);
	//Swaption swaption1(swap, true);
	//Swaption swaption2(swap, false);


	//std::cout << "=== Formule fermee ==="  << "\n";
	//std::cout << "Payeur \t\t: " << percentage(pricer_swap.Evaluate(swap)) << "\n";
	//std::cout << "Payeuse \t: " << percentage(pricer_option1.Evaluate(swaption1)) << "\n";
	//std::cout << "Receveuse \t: " << percentage(pricer_option1.Evaluate(swaption2)) << "\n";
	//std::cout << "=== Arbre trinomial ==="  << "\n";
	//std::cout << "Payeuse \t: " << percentage(pricer_option2.Evaluate(swaption1)) << "\n";
	//std::cout << "Receveuse \t: " << percentage(pricer_option2.Evaluate(swaption2)) << "\n";

	//std::cout << "\n" << pricer_option2;
	std::ofstream(std::ofstream("tree.gv", std::ios::trunc)) << pricer_option2;


	////std::ofstream(std::ofstream("arbre.csv", std::ios::trunc)) << arbre;

	std::getchar();
	return 0;
}