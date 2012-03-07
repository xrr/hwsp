#pragma once
#include <iostream>
#include <string>
#include <vector>

#include <math.h>
#include <boost/math/special_functions/round.hpp>

#include "Tools.h"


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

class RateCurve {
public:
	std::vector<double> times;
	std::vector<double> rates;
	double zerocoupon(double time) {
		double prev_time = 0, tmp = ((time<=times[0])?rates[0]*time:0);
		unsigned int k=0; while((k<times.size()) && (times[k]<time)) {
			tmp+=rates[k]*(std::min<double>(time,((k==times.size()-1)?time:times[k+1]))-prev_time);
			if (times.size()-1!=k) prev_time = times[k+1];
			k++;
		};
		return exp(-tmp);
	};
	double rate(double time) {
		if (time<=times[0]) return rates[0]; 
		else {
			unsigned int k=0; while((k<times.size()) && (times[k]<=time)) k++;
			return rates[k-1];
		}
	};
};

class HullWhite {
public:
	double a, sigma;
	HullWhite(double _a, double _sigma) : a(_a), sigma(_sigma) {};
};

class PricerGeneric {
public :
	RateCurve ratecurve;
	PricerGeneric(RateCurve _ratecurve) : ratecurve(_ratecurve) {};

	virtual double Evaluate(Swap swap) {
		double PrixSwap=0;
		double prev_paymentdate = swap.startdate;
		for (unsigned int i=0; i<=swap.paymentdates.size()-1; ++i) {
			PrixSwap += ratecurve.zerocoupon(swap.paymentdates[i])*( // Discount Factor
				swap.strikerate*(swap.paymentdates[i]-prev_paymentdate) // Fixed Leg
				- (ratecurve.zerocoupon(prev_paymentdate)/ratecurve.zerocoupon(swap.paymentdates[i])-1)); // Float Leg
			prev_paymentdate = swap.paymentdates[i];
		}
		return PrixSwap;
	};
};

class PricerOptions : public PricerGeneric {
public:
	HullWhite hullwhite;
	PricerOptions(RateCurve _ratecurve, HullWhite _hullwhite) : PricerGeneric(_ratecurve), hullwhite(_hullwhite) {};
	virtual double Evaluate(Swaption)=0;
};

class ClosedFormula;

class ClosedFormulaSwaptionFunctor {
private :
	ClosedFormula *pCF_;
	Swaption swaption_;
public :
	ClosedFormulaSwaptionFunctor(ClosedFormula *, Swaption);
	double operator()(double);
};

class ClosedFormula : public PricerOptions {

public :

	double B(double t, double T) {return 1/hullwhite.a*(1-exp(-hullwhite.a*(T-t)));};

	double A(double t, double T) {
		return ratecurve.zerocoupon(T)/ratecurve.zerocoupon(t)*exp(
			B(t,T)*ratecurve.rate(t)-hullwhite.sigma*hullwhite.sigma/(4*hullwhite.a)*(1-exp(-2*hullwhite.a)*B(t,T)*B(t,T)));
	};
	ClosedFormula(RateCurve _ratecurve, HullWhite _hullwhite) : PricerOptions(_ratecurve, _hullwhite) {};

	double f(Swaption swaption, double r) {
		double res=0;
		double prev_paymentdate = swaption.swap.startdate;
		for(std::vector<double>::size_type i=0; i<=swaption.swap.paymentdates.size()-1; ++i) { // todo
			double c=swaption.swap.strikerate*(swaption.swap.paymentdates[i]-prev_paymentdate);
			prev_paymentdate = swaption.swap.paymentdates[i];
			if (i==swaption.swap.paymentdates.size()-1) ++c;
			double X = A(swaption.swap.startdate, swaption.swap.paymentdates[i])*exp(-B(swaption.swap.startdate,swaption.swap.paymentdates[i])*r);
			res += c*X;
		}
		return res;
	};

	virtual double Evaluate(Swaption swaption) {

		double r_star = dichotomicRootSearch<double, double, ClosedFormulaSwaptionFunctor>
			(ClosedFormulaSwaptionFunctor(this,swaption), 1.0, 0.00001, -1.00, 1.00);

		int w = (swaption.ispayeroption?1:-1);
		double PrixSwaption = 0;
		double prev_paymentdate = swaption.swap.startdate;
		for(std::vector<double>::size_type i=0; i<=swaption.swap.paymentdates.size()-1; ++i) {
			double c = swaption.swap.strikerate*(swaption.swap.paymentdates[i]-prev_paymentdate);
			if (swaption.swap.paymentdates.size()-1==i) ++c;
			prev_paymentdate = swaption.swap.paymentdates[i];
			double X = A(swaption.swap.startdate, swaption.swap.paymentdates[i])*exp(-B(swaption.swap.startdate,swaption.swap.paymentdates[i])*r_star);
			double sigma_p = hullwhite.sigma
				* B(swaption.swap.startdate,swaption.swap.paymentdates[i])
				* sqrt((1-exp(-2*hullwhite.a*(swaption.swap.startdate)))/(2*hullwhite.a));
			double h = sigma_p/2 + (1/sigma_p) * log(
				ratecurve.zerocoupon(swaption.swap.paymentdates[i])
				/(ratecurve.zerocoupon(swaption.swap.startdate)*X));

			BoostGauss gauss;
			double ZBO = w * X * ratecurve.zerocoupon(swaption.swap.startdate) * gauss.cdf(w*(sigma_p-h))
				- w * ratecurve.zerocoupon(swaption.swap.paymentdates[i]) * gauss.cdf(-w*h);
			PrixSwaption += c*ZBO;
		}
		return PrixSwaption;
	};

};

ClosedFormulaSwaptionFunctor::ClosedFormulaSwaptionFunctor(ClosedFormula *_pCF, Swaption _swaption) : pCF_(_pCF), swaption_(_swaption) {};
double ClosedFormulaSwaptionFunctor::operator()(double x) {return pCF_->f(swaption_,x);};

class Node;

class Transition {
public:
	double probability;
	Node *destination;
	Transition(void) : probability(NULL), destination(NULL) {};
	Transition(double _probability) : probability(_probability), destination(NULL) {};
};

class Node {
public:
	int relative_position;
	double x, alpha, q; // driftless rate, displacement, arrow_debreu
	double value; // todo : CHANGE
	Transition transitions[3]; // down, middle, up
	double r(void) { return x+alpha; }; // displaced rate
	Node(int _relative_position, double _x) : relative_position(_relative_position), x(_x), alpha(NULL), q(0), value(0) {};
	Node(int _relative_position, double _x, double _q, double _alpha) : relative_position(_relative_position), x(_x), alpha(_alpha), q(_q), value(0) {};
};



// todo = tree destructor
// todo = clean X= twice
class Tree : public PricerOptions {

private:
	std::vector<std::vector<Node*>> slices_;
	std::vector<double> dates_;




public:

	void resetCashflowValues(void) {
		for (auto pslice = slices_.begin(); pslice < slices_.end(); ++pslice)
			for (auto ppnode = pslice->begin(); ppnode < pslice->end(); ++ppnode)
				(*ppnode)->value = 0;
	};

	friend std::ostream& operator<< (std::ostream&, Tree&);

	virtual double Evaluate(Swaption swaption) {
		//todo : check dates if dates of swap = last dates of tree
		resetCashflowValues();
		auto pdate = dates_.rbegin(); auto pslice = slices_.rbegin();
		while (pdate < dates_.rend()) {
			// DISCLAIMER /!\ *(pdate+1) = t-1 : REVERSE iterators
			for (auto ppnode = pslice->begin(); ppnode < pslice->end(); ++ppnode) {
				if (pdate != dates_.rbegin()) for (int l=0; l<3; ++l) // discounted & probabilized value of each successor (except for last slice in tree)
					((*ppnode)->value) += exp(-(*ppnode)->r()*(*(pdate-1)-*pdate))*(*ppnode)->transitions[l].probability*(*ppnode)->transitions[l].destination->value;
				if (*pdate>swaption.swap.startdate)
					((*ppnode)->value) += swaption.swap.strikerate*(*pdate-*(pdate+1))-(exp((*ppnode)->r()*(*pdate-*(pdate+1)))-1); // fixed - float 
				if (*pdate == swaption.swap.startdate)
					((*ppnode)->value) = std::max<double>(0,(swaption.ispayeroption?-1:1)*(*ppnode)->value); // option exercice todo : max(0,+/- cashflow if receiver swaption)
				}
			++pdate; ++pslice;
		}
		for (auto ppnode = pslice->begin(); ppnode < pslice->end(); ++ppnode)
				if (pdate != dates_.rbegin()) for (int l=0; l<3; ++l) // last computation to go back to time zero
					((*ppnode)->value) += exp(- (*ppnode)->r() * (*(pdate-1)-0)) * (*ppnode)->transitions[l].probability * (*ppnode)->transitions[l].destination->value;
		return (*pslice)[0]->value; // NPV
	};

	//todo : error returned if dates are not in right order, and if first date is <=0 (=0 included, it will fail)
	Tree (RateCurve _ratecurve, HullWhite _hullwhite, std::vector<double> _dates) : PricerOptions(_ratecurve, _hullwhite), dates_(_dates) {
		slices_.push_back(std::vector<Node*>(1,new Node(0,0,1,ratecurve.rates[0]))); // todo : remplacer rates[0] vraiment bizarre
		double prev_date = 0;
		for (auto pdate = dates_.begin(); pdate < dates_.end(); ++pdate) {
			slices_.push_back(std::vector<Node*>());
			auto pnewslice = slices_.rbegin(); // last slice
			auto pcurrentslice = &(*(pnewslice+1)); // previous slice
			double delta_t = *pdate - prev_date;
			double e = exp(-hullwhite.a*delta_t);
			double V = hullwhite.sigma * sqrt((1-e*e)/(2*hullwhite.a));
			double delta_x = V * sqrt(3.0);
			for (auto ppcurrentnode = pcurrentslice->begin(); ppcurrentnode < pcurrentslice->end(); ++ppcurrentnode) {
				double M = (*ppcurrentnode)->x*e;			
				int k = boost::math::iround(M/delta_x);
				double R = (M-k*delta_x)/V; // eta/V
				for (int l=-1; l<=+1; ++l) {
					Transition newtransition((1+R*R)/6.0+((0==l)?(1-R*R):l*R/sqrt(3.0))/2.0);
					auto ppexistingnode = pnewslice->begin();
					while (ppexistingnode < pnewslice->end() && newtransition.destination == NULL) {
						if ((*ppexistingnode)->relative_position == k+l)
							newtransition.destination = *ppexistingnode;
						++ppexistingnode;
					}
					if (newtransition.destination == NULL) {
						pnewslice->push_back(new Node(k+l, (k+l)*delta_x));
						newtransition.destination = *(pnewslice->rbegin());
					}
					(*ppcurrentnode)->transitions[1+l]=newtransition;
					newtransition.destination->q+=(*ppcurrentnode)->q*newtransition.probability*exp(-(*ppcurrentnode)->r()*delta_t);
				}
			}
			double alpha = 0; // displacement
			for (auto ppcurrentnode = pnewslice->begin(); ppcurrentnode != pnewslice->end(); ++ppcurrentnode)
				alpha+=(*ppcurrentnode)->q*exp(-(*ppcurrentnode)->x*delta_t);
			double tplus1 = (pdate+1<dates_.end())?*(pdate+1):*pdate+delta_t;
			alpha = log(alpha/(ratecurve.zerocoupon(tplus1)))/(tplus1-*pdate);
			for (auto ppcurrentnode = pnewslice->begin(); ppcurrentnode != pnewslice->end(); ++ppcurrentnode)
				(*ppcurrentnode)->alpha = alpha;
			prev_date = *pdate;
		}
	};
};

//Graphviz output (see http://en.wikipedia.org/wiki/DOT_language for eg.)
std::ostream& operator<< (std::ostream& flux, Tree& tree) {
	flux << "digraph Tree {" << std::endl
		<< "graph [rankdir=\"LR\",splines=false,label=\"Trinomial Tree\"];" << std::endl
		<< "node [shape=record,color=blue];" << std::endl
		<< "edge [style=dashed,color=red];" << std::endl
		<< std::endl;
	for (auto pcurrentslice = tree.slices_.begin(); pcurrentslice < tree.slices_.end(); ++pcurrentslice) {
		for (auto ppcurrentnode = pcurrentslice->rbegin(); ppcurrentnode < pcurrentslice->rend(); ++ppcurrentnode) {
			flux <<  "n" << (*ppcurrentnode) << " [label=\""
				<< "short rate : " << percentage((*ppcurrentnode)->r()) << "\\n "
				<< "driftless rate : " << percentage((*ppcurrentnode)->x) << "\\n "
				<< "contingent claim PV : " << percentage((*ppcurrentnode)->q) << "\\n\\n "
				<< "cashflow value : " << percentage((*ppcurrentnode)->value)
				<< "\"];"<< std::endl;
			for (int l=2; 0<=l; --l)
				if ((*ppcurrentnode)->transitions[l].destination != NULL)
					flux << "n" << (*ppcurrentnode) << " -> n" << (*ppcurrentnode)->transitions[l].destination
					<< " [label=\"P(" << (0<l?(1<l?"up":"mid."):"down") << ") = " << percentage((*ppcurrentnode)->transitions[l].probability) << "\"];"<< std::endl;
		} flux << std::endl;
	} flux << "}" << std::endl;
	return flux;
};