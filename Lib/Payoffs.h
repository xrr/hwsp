#pragma once
#include <vector>

class Dates {
public:
	std::vector<double> vect;
	//private :
	//	std::vector<double> vect_;
	//public :
	Dates(std::vector<double> values) : vect(values) {};
	double operator()(int i) {return vect[i];};
	int size(void) {return vect.size();};
};

class Payoff {};

// by convention : payer swap (receives float and pays fixed in exchange)
class Swap : public Payoff {
public :
	double strikerate;
	double startdate;
	Dates paymentdates;
	Swap(double _strikerate, double _startdate, Dates _paymentdates) : strikerate(_strikerate), startdate(_startdate), paymentdates(_paymentdates) {};
};

class Swaption : public Payoff {
public :
	Swap swap;
	bool ispayeroption; // false : receiver swaption
	Swaption(Swap _swap, Dates _exercisedates, bool _ispayeroption) :
		swap(_swap), ispayeroption(_ispayeroption) {};
};