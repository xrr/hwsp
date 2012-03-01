#include <iostream>
#include "Node.h"

#include <math.h>

int main() {

	Curve courbe;

	for (int i = 0 ; i<3 ; i++) courbe.times.push_back(1+i);
	courbe.rates.push_back(0.0382365);
	courbe.rates.push_back(0.0451162);
	courbe.rates.push_back(0.0508626);

	//for (int i = 0 ; i<4 ; i++) courbe.times.push_back(1+i);
	//courbe.rates.push_back(0.05092755);
	//courbe.rates.push_back(0.05795397);
	//courbe.rates.push_back(0.06304557);
	//courbe.rates.push_back(0.06733466);

	HullWhite hw;
	hw.a = 0.1;
	hw.sigma = 0.01;
	hw.ratecurve = courbe;

	Tree arbre(hw);
	arbre.construct();
	std::cout << "\n" << arbre;
	std::getchar();
}