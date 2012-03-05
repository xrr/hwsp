#pragma once
#include <vector>


class Payoff {};

// by convention : payer swap (receives float and pays fixed in exchange)
class Swap : public Payoff {
public :
	double strikerate;
	double startdate;
	std::vector<double> paymentdates;
	Swap(double _strikerate, double _startdate, std::vector<double> _paymentdates) : strikerate(_strikerate), startdate(_startdate), paymentdates(_paymentdates) {};
};

class Swaption : public Payoff {
public :
	Swap swap;
	bool ispayeroption; // false : receiver swaption
	Swaption(Swap _swap, bool _ispayeroption) :
		swap(_swap), ispayeroption(_ispayeroption) {};
};