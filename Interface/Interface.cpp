#include <iostream>
#include <cmath>
#include "Gauss.h"

int main() {
	Gauss* myCDF1 = new AbramowitzStegunGauss;
	Gauss* myCDF2 = new BoostGauss;
	std::cout << myCDF1->erf(-1) << " ; " << myCDF1->cdf(1) << std::endl;
	std::cout << myCDF2->erf(-1) << " ; " << myCDF2->cdf(1) << std::endl;
	//std::getchar();
	delete myCDF1, myCDF2;
}