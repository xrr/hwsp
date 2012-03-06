#pragma once
#include <iostream>
#include <string>
#include <vector>

#include <math.h>
#include <boost/math/special_functions/round.hpp>

#include "Models.h"
#include "Payoffs.h"
#include "Tools.h"
#include "Gauss.h"

//#include <algorithm>

//todo = check transition pointers in class Tree

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
		for(std::vector<double>::size_type i=0; i<=swaption.swap.paymentdates.size()-1; ++i) { // todo
			double prev_paymentdate = swaption.swap.startdate;
			double c=swaption.swap.strikerate*(swaption.swap.paymentdates[i]-prev_paymentdate);
			prev_paymentdate = swaption.swap.paymentdates[i];
			if (i=swaption.swap.paymentdates.size()-1) {c=c+1;}
			double X = A(swaption.swap.startdate, swaption.swap.paymentdates[i])*exp(-B(swaption.swap.startdate,swaption.swap.paymentdates[i])*r);
			res += c*X;
		}
		return res;
	};

	virtual double Evaluate(Swaption swaption) {

		//dichotomic root search of r*
		double f_target = 1;
		double min=-1, max=1, target;
		double f_min = f(swaption, min)-f_target, f_max = f(swaption, max)-f_target;
		double precision = 0.00001, f_target_temp = abs(precision)+1234; // arbitrary >0 constant to enter into while loop
		while (abs(f_target_temp) > precision) { // seeking for function's zero
			target = (min+max) / 2 ;
			f_target_temp = f(swaption, target)-f_target;
			if (f_min*f_target_temp < 0) { 
				max = target;
				f_max = f_target_temp;
			} else  {
				min = target;
				f_max = f_target_temp;
			}
		}

		double r_star = target;

		double PrixSwaption = 0;
		for(std::vector<double>::size_type i=0; i<=swaption.swap.paymentdates.size()-1; ++i) {
			double prev_paymentdate = swaption.swap.startdate;
			double c = swaption.swap.strikerate*(swaption.swap.paymentdates[i]-prev_paymentdate);
			if (swaption.swap.paymentdates.size()-1==i) ++c; // ??
			prev_paymentdate = swaption.swap.paymentdates[i]; // ??
			double X = A(swaption.swap.startdate, swaption.swap.paymentdates[i])*exp(-B(swaption.swap.startdate,swaption.swap.paymentdates[i])*r_star);
			double sigma_p = hullwhite.sigma
				* B(swaption.swap.startdate,swaption.swap.paymentdates[i])
				* sqrt((1-exp(-2*hullwhite.a*(swaption.swap.startdate)))/(2*hullwhite.a));
			double h = sigma_p/2 + (1/sigma_p) * log(
				ratecurve.zerocoupon(swaption.swap.paymentdates[i])
				/(ratecurve.zerocoupon(swaption.swap.startdate)*X));
			int w = (swaption.ispayeroption?1:-1);
			AbramowitzStegunGauss gauss;
			double ZBO = w * X * ratecurve.zerocoupon(swaption.swap.startdate) * gauss.cdf(w*(sigma_p-h))
				- w * ratecurve.zerocoupon(swaption.swap.paymentdates[i]) * gauss.cdf(-w*h);
			PrixSwaption += c*ZBO;
		}
		return PrixSwaption;
	};

// Romain test, with an external rootsearch method
//	virtual double Evaluate2(Swaption swaption) {
//		inline g(double) ;
//		g(r) = f(swaption,r_);
//		typedef double mfunc(Swaption,double);
//		double r_star = dichotomicRootSearch<double, double, mfunc>(&ClosedFormula::f, 1, 0.0010, 0, 1);
//		//double r_star = dichotomicRootSearch<double, double>(1, 0.0010, 0,1);
//
//		double PrixSwaption = 0;
//		for(std::vector<double>::size_type i=0; i<=swaption.swap.paymentdates.size()-1; ++i) {
//			double prev_paymentdate = swaption.swap.startdate;
//			double c = swaption.swap.strikerate*(swaption.swap.paymentdates[i]-prev_paymentdate);
//			if (swaption.swap.paymentdates.size()-1==i) ++c; // ??
//			prev_paymentdate = swaption.swap.paymentdates[i]; // ??
//			double X = A(swaption.swap.startdate, swaption.swap.paymentdates[i])*exp(-B(swaption.swap.startdate,swaption.swap.paymentdates[i])*r_star);
//			double sigma_p = hullwhite.sigma
//				* B(swaption.swap.startdate,swaption.swap.paymentdates[i])
//				* sqrt((1-exp(-2*hullwhite.a*(swaption.swap.startdate)))/(2*hullwhite.a));
//			double h = sigma_p/2 + (1/sigma_p) * log(
//				ratecurve.zerocoupon(swaption.swap.paymentdates[i])
//				/(ratecurve.zerocoupon(swaption.swap.startdate)*X));
//			int w = (swaption.ispayeroption?1:-1);
//			AbramowitzStegunGauss gauss;
//			double ZBO = w * X * ratecurve.zerocoupon(swaption.swap.startdate) * gauss.cdf(w*(sigma_p-h))
//				- w * ratecurve.zerocoupon(swaption.swap.paymentdates[i]) * gauss.cdf(-w*h);
//			PrixSwaption += c*ZBO;
//		}
//		return PrixSwaption;
//	};
};

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
	double x, alpha, q; // no drift rate, displacement, arrow_debreu
	double value; // todo : CHANGE
	Transition transitions[3]; // down, middle, up
	double r(void) { return x+alpha; }; // displaced rate
	Node(int _relative_position, double _x) : relative_position(_relative_position), x(_x), alpha(NULL), q(0), value(0) {};
	Node(int _relative_position, double _x, double _q, double _alpha) : relative_position(_relative_position), x(_x), alpha(_alpha), q(_q), value(0) {};
};

class Tree : public PricerOptions {

private:
	std::vector<std::vector<Node>> slices_;
	std::vector<double> dates_;

	void resetCashflowValues(void) {
		for (auto pslice = slices_.begin(); pslice != slices_.end(); ++pslice)
			for (auto pnode = pslice->begin(); pnode != pslice->end(); ++pnode)
				pnode->value = 0;
	};

public:

	friend std::ostream& operator<< (std::ostream&, Tree&);

	virtual double Evaluate(Swaption swaption) {
		//todo : check dates if dates of swap = last dates of tree
		resetCashflowValues();
		auto pdate = dates_.end(); auto pslice = slices_.end();
		while (pdate >= dates_.begin()) {
			for (auto pnode = pslice->begin(); pnode != pslice->end(); ++pnode) {
				if (*pdate>swaption.swap.startdate)
					pnode->value += swaption.swap.strikerate * (*pdate-*(pdate-1)) - (exp(pnode->r()*(*pdate-*(pdate-1)))-1); //fixed - float 
				if (pdate != ratecurve.times.end()) for (int l=0; l<3; ++l) // discounted & probabilized value of each successor (except for last slice in tree)
					pnode->value += exp(- pnode->r() * (*(pdate+1)-*pdate)) * pnode->transitions[l].probability * pnode->transitions[l].destination->value; 
				if (*pdate = swaption.swap.startdate)
					pnode->value = std::max<double>(0,(swaption.ispayeroption?-1:1)*pnode->value); // option exercice todo : max(0,+/- cashflow if receiver swaption)
			}
			--pdate; --pslice;
		}
		return (*pslice)[0].value; // NPV
	};

	Tree (RateCurve _ratecurve, HullWhite _hullwhite, std::vector<double> dates) : PricerOptions(_ratecurve, _hullwhite) {
		slices_.push_back(std::vector<Node>(1,Node(0,0,1,ratecurve.rates[0])));
		//std::vector<Node> slice0;
		//slice0.push_back(Node(0,0,1,ratecurve.rates[0]));
		//slices.push_back(slice0);
		double prev_date = 0;
		for (auto pdate = dates.begin(); pdate != dates.end()-1; ++pdate) {
			std::vector<Node> newslice;
			double delta_t = *pdate - prev_date;
			double e = exp(-hullwhite.a*delta_t);
			double V = hullwhite.sigma * sqrt((1-e*e)/(2*hullwhite.a));
			double delta_x = V * sqrt(3.0);
			for (auto currentnode = slices_.back().begin(); currentnode != slices_.back().end(); ++currentnode) {
				double M = currentnode->x*e;			
				int k = boost::math::iround(M/delta_x);
				double R = (M-k*delta_x)/V; // eta/V
				for (int l=-1; l<=+1; ++l) {
					Transition newtransition((1+R*R)/6.0+((0==l)?(1-R*R):l*R/sqrt(3.0))/2.0);
					auto existingnode = newslice.begin();
					while (existingnode != newslice.end() && newtransition.destination == NULL) {
						if (existingnode->relative_position == k+l) newtransition.destination = &(*existingnode);
						++existingnode;
					}
					if (newtransition.destination == NULL) {
						newslice.push_back(Node(k+l, (k+l)*delta_x));
						newtransition.destination=&newslice.back();
					}
					currentnode->transitions[1+l]=newtransition;
					newtransition.destination->q+=currentnode->q*newtransition.probability*exp(-currentnode->r()*delta_t);
				}
			}
			slices_.push_back(newslice);
			double alpha = 0;
			for (auto currentnode = slices_.back().begin(); currentnode != slices_.back().end(); ++currentnode)
				alpha+=currentnode->q*exp(-currentnode->x*delta_t);
			alpha = log(alpha/(ratecurve.zerocoupon(*(pdate+1))))/delta_t;
			std::cout << fixed_percentage(alpha) << "\n";
			for (auto currentnode = slices_.back().begin(); currentnode != slices_.back().end(); ++currentnode)
				currentnode->alpha = alpha;
			prev_date = *pdate;
		}
	};

};

std::ostream& operator<< (std::ostream& flux, Tree& tree) {
	for (auto currentslice = tree.slices_.begin(); currentslice != tree.slices_.end(); ++currentslice) {
		for (auto currentnode = currentslice->begin(); currentnode != currentslice->end(); ++currentnode) {
			flux << fixed_percentage(currentnode->r(),5) << "; ";
			/*flux << fixed_percentage(currentnode->transitions[2].probability) << " "
			<< fixed_percentage(currentnode->transitions[1].probability) << " "
			<< fixed_percentage(currentnode->transitions[0].probability) << " ";
			flux << "\n";*/
		}
		flux << "\n";
	}
	return flux;
}