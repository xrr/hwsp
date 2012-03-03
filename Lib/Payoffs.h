#pragma once
#include <vector>

class Dates {
public:
	std::vector<double> vect;
	//private :
	//	std::vector<double> vect_;
	//public :
	//	Dates(std::vector<double> values) : vect_(values) {};
	double operator()(int i) {return vect[i];};
	int size(void) {return vect.size();};
};

class Payoff {};

// by convention : payer swap (receives float and pays fixed in exchange)
class Swap : public Payoff {
public :
	double strikerate;
	Dates paymentdates;
	Swap(double _strikerate, Dates _paymentdates) : strikerate(_strikerate), paymentdates(_paymentdates) {};
};

class Swaption : public Payoff {
public :
	Swap swap;
	Dates exercisedates; // first date implementation for now (no berm option)
	bool ispayeroption; // false : receiver swaption
	Swaption(Swap _swap, Dates _exercisedates, bool _ispayeroption) :
		swap(_swap), exercisedates(_exercisedates), ispayeroption(_ispayeroption) {};
};