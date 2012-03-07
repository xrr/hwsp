#include <iostream>
#include <fstream>
#include <vector>

#include <math.h>

#include "Lib.h"



int main () {

	RateCurve courbe;
	for (int i = 0 ; i<3 ; i++) {
		courbe.times.push_back(2+i);
		courbe.rates.push_back(0.05);
	}

	std::vector<double> tree_dates;
	for (int i = 0 ; i<5 ; i++) tree_dates.push_back(0.2*(1+i));
	tree_dates.insert(tree_dates.end(), courbe.times.begin(), courbe.times.end());

	HullWhite hullwhite(0.1, 0.01);


	PricerGeneric pricer_swap(courbe);
	ClosedFormula pricer_option1(courbe, hullwhite);
	Tree pricer_option2(courbe,hullwhite, tree_dates);

	for (int i = 0; i < 5; ++i) {

		double strike = 0.03+0.01*i;

		Swap swap(strike,1,courbe.times);
		Swaption swaption1(swap, true);
		Swaption swaption2(swap, false);

		std::cout
			<< std::endl << "=== Payer Swap ===" 
			<< std::endl << "Strike \t: " << percentage(strike)
			<< std::endl << "Price  \t: " << percentage(pricer_swap.Evaluate(swap))
			<< std::endl
			<< std::endl << "=== Payer Option === "
			<< std::endl << "Formula\t: " << percentage(pricer_option1.Evaluate(swaption1))
			<< std::endl << "Tree   \t: " << percentage(pricer_option2.Evaluate(swaption1))
			<< std::endl
			<< std::endl << "=== Rec. Option === "
			<< std::endl << "Formula\t: " << percentage(pricer_option1.Evaluate(swaption2))
			<< std::endl << "Tree   \t: " << percentage(pricer_option2.Evaluate(swaption2))
			<< std::endl << std::endl << std::endl << std::endl << std::endl;
	}

	std::ofstream(std::ofstream("tree.gv", std::ios::trunc)) << pricer_option2;

	std::getchar();
	return 0;
}