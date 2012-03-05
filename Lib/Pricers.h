#pragma once
#include <iostream>
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
		for (int i=0; i<=swap.paymentdates.size()-1; i++) {
			PrixSwap += swap.strikerate*ratecurve.zerocoupon(swap.paymentdates(i)) - 
				(1/(swap.paymentdates(i)- prev_paymentdate)*(ratecurve.zerocoupon(i)/ratecurve.zerocoupon(prev_paymentdate)-1))*ratecurve.zerocoupon(swap.paymentdates(i)) ;
			prev_paymentdate = swap.paymentdates(i);
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
	Transition transitions[3]; // down, middle, up
	double r(void) { return x+alpha; }; // displaced rate
	Node(int _relative_position, double _x) : relative_position(_relative_position), x(_x), alpha(NULL), q(0) {};
	Node(int _relative_position, double _x, double _q, double _alpha) : relative_position(_relative_position), x(_x), alpha(_alpha), q(_q) {};
};

class Tree : public PricerOptions {

private:
	std::vector<std::vector<Node>> slices;

public:

	friend std::ostream& operator<< (std::ostream&, Tree&);

	virtual double Evaluate(Swaption swaption) {
		//todo;
		return -999;
	};


	Tree (RateCurve _ratecurve, HullWhite _hullwhite) : PricerOptions(_ratecurve, _hullwhite) {
		std::vector<Node> slice0;
		slice0.push_back(Node(0,0,1,ratecurve.rates[0]));
		slices.push_back(slice0);
		double prev_time = 0;
		for (auto time = ratecurve.times.begin(); time != ratecurve.times.end()-1; time++) {
			std::vector<Node> newslice;
			double delta_t = *time-prev_time;
			double e = exp(-hullwhite.a*delta_t);
			double V = hullwhite.sigma * sqrt((1-e*e)/(2*hullwhite.a));
			double delta_x = V * sqrt(3.0);
			for (auto currentnode = slices.back().begin(); currentnode != slices.back().end(); currentnode++) {
				double M = currentnode->x*e;			
				int k = boost::math::iround(M/delta_x);
				double R = (M-k*delta_x)/V; // eta/V
				for (int l=-1; l<=+1; l++) {
					Transition newtransition((1+R*R)/6.0+((0==l)?(1-R*R):l*R/sqrt(3.0))/2.0);
					auto existingnode = newslice.begin();
					while (existingnode != newslice.end() && newtransition.destination == NULL) {
						if (existingnode->relative_position == k+l) newtransition.destination = &(*existingnode);
						existingnode++;
					}
					if (newtransition.destination == NULL) {
						newslice.push_back(Node(k+l, (k+l)*delta_x));
						newtransition.destination=&newslice.back();
					}
					currentnode->transitions[1+l]=newtransition;
					newtransition.destination->q+=currentnode->q*newtransition.probability*exp(-currentnode->r()*delta_t);
				}
			}
			slices.push_back(newslice);
			double alpha = 0;
			for (auto currentnode = slices.back().begin(); currentnode != slices.back().end(); currentnode++) alpha+=currentnode->q*exp(-currentnode->x*delta_t);
			alpha = log(alpha/(ratecurve.zerocoupon(*(time+1))))/delta_t;
			std::cout << fixed_percentage(alpha) << "\n";
			for (auto currentnode = slices.back().begin(); currentnode != slices.back().end(); currentnode++) currentnode->alpha = alpha;
			prev_time = *time;
		}
	};

};

std::ostream& operator<< (std::ostream& flux, Tree& tree) {
	for (auto currentslice = tree.slices.begin(); currentslice != tree.slices.end(); currentslice++) {
		for (auto currentnode = currentslice->begin(); currentnode != currentslice->end(); currentnode++) {
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


class ClosedFormula : public PricerOptions {
private :
	double B(double t, double T) {return 1/hullwhite.a*(1-exp(-hullwhite.a*(T-t)));};
	double A(double t, double T) {
		return ratecurve.zerocoupon(T)/ratecurve.zerocoupon(t)*exp(
			B(t,T)*ratecurve.rate(t)-hullwhite.sigma*hullwhite.sigma/(4*hullwhite.a)*(1-exp(-2*hullwhite.a)*B(t,T)*B(t,T)));
	};
	double f(Swaption swaption, double r) {
		double res=0;
		for(int i=0; i<=swaption.swap.paymentdates.size()-1; i++) {
					double prev_paymentdate = swaption.swap.startdate;
					double c=swaption.swap.strikerate*(swaption.swap.paymentdates(i)-prev_paymentdate);
					prev_paymentdate = swaption.swap.paymentdates(i);
					if (i=swaption.swap.paymentdates.size()-1) {c=c+1;}
					double X = A(swaption.swap.startdate, swaption.swap.paymentdates(i))*exp(-B(swaption.swap.startdate,swaption.swap.paymentdates(i))*r);
				res += c*X;
				}
	return res;
	};

public :
	ClosedFormula(RateCurve _ratecurve, HullWhite _hullwhite) : PricerOptions(_ratecurve, _hullwhite) {};

	virtual double Evaluate(Swaption swaption) {
		//Trouver r_star par dichotomie
		
		double f_target = 1, precision = 0.00001;
		double f_target_temp = abs(precision)+1234; // 1234 : arbitrary positive constant to enter into while loop

		double min=-1, max=1, target;

		double f_min = f(swaption, min)-f_target;
		double f_max = f(swaption, max)-f_target;

		while (abs(f_target_temp) > precision) { // seeking for function's zero

			target = (min+max) / 2 ;
			f_target_temp = f(swaption, target)-f_target;
			
			if (f_min*f_target_temp < 0) { 
				max = target;
				f_max = f_target_temp;
			}
			else  {
				min = target;
				f_max = f_target_temp;
			}

		}
		
		double r_star = target;
		//std::cout << r_star << " -> " << f(swaption,r_star) << "\n";
	
		
		AbramowitzStegunGauss gauss;
		double PrixSwaption_payer = 0, PrixSwaption_receiver=0;
		for(int i=0; i<=swaption.swap.paymentdates.size()-1; i++) {
			double prev_paymentdate = swaption.swap.startdate; // remplacer 0 par startdate
			double c=swaption.swap.strikerate*(swaption.swap.paymentdates(i)-prev_paymentdate);
			prev_paymentdate = swaption.swap.paymentdates(i);
			if (i=swaption.swap.paymentdates.size()-1) {c=c+1;}
			double X = A(swaption.swap.startdate, swaption.swap.paymentdates(i))*exp(-B(swaption.swap.startdate,swaption.swap.paymentdates(i))*r_star);
			double sigma_p = hullwhite.sigma*sqrt((1-exp(-2*hullwhite.a*(swaption.swap.startdate-swaption.swap.paymentdates(i))))/(2*hullwhite.a))*B(swaption.swap.startdate,swaption.swap.paymentdates(i));
			double h = (1.0/sigma_p)*log(ratecurve.zerocoupon(swaption.swap.paymentdates(i))/(ratecurve.zerocoupon(swaption.swap.startdate)*X))+sigma_p/2;
			double ZBP = X*ratecurve.zerocoupon(swaption.swap.startdate)*gauss.cdf(-h+sigma_p) - ratecurve.zerocoupon(swaption.swap.paymentdates(i))*gauss.cdf(-h);
			double ZBC = -X*ratecurve.zerocoupon(swaption.swap.startdate) + ratecurve.zerocoupon(swaption.swap.paymentdates(i)) +ZBP;
			PrixSwaption_payer += c*ZBP;
			PrixSwaption_receiver += c*ZBC;
		}//boucle swaption
		if (swaption.ispayeroption == true) {return PrixSwaption_payer;}
		else {return PrixSwaption_receiver;}
	};
};


