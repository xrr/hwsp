#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <boost/math/special_functions/round.hpp>
#include <sstream>
#include <algorithm>

//todo = check transition pointers in class Tree

std::string fixed_percentage (double value, int precision = 3) {
	std::ostringstream out (std::ostringstream::out);
	if (0==value) precision--;
	out.precision(precision);
	out << std::showpoint << std::showpos << std::setw(precision) << std::left;
	out << 100*value << "%";
	return out.str();
};

class Curve {
public:
	std::vector<double> times;
	std::vector<double> rates;
	double zerocoupon(double time) {
		double prev_time = 0, tmp = 0;
		unsigned int k=0; while((k < times.size()) && (times[k]<=time)) {
			tmp+=rates[k]*(times[k]-prev_time);
			prev_time = times[k];
			k++;
		};
		return exp(-tmp);
	};
};

class HullWhite {
public:
	double a, sigma;
	Curve ratecurve;
	HullWhite(double _a, double _sigma, Curve _ratecurve) : a(_a), sigma(_sigma), ratecurve(_ratecurve) {};
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

class Tree {

private:
	HullWhite hw_;
	std::vector<std::vector<Node>> dates;

public:

	friend std::ostream& operator<< (std::ostream&, Tree&);

	Tree (HullWhite hw) : hw_(hw) {
		std::vector<Node> date0;
		date0.push_back(Node(0,0,1,hw_.ratecurve.rates[0]));
		dates.push_back(date0);
	};

	void construct(void) {
		double prev_time = 0;
		for (auto time = hw_.ratecurve.times.begin(); time != hw_.ratecurve.times.end()-1; time++) {
			std::vector<Node> newdate;
			double delta_t = *time-prev_time;
			double e = exp(-hw_.a*delta_t);
			double V = hw_.sigma * sqrt((1-e*e)/(2*hw_.a));
			double delta_x = V * sqrt(3.0);
			for (auto currentnode = dates.back().begin(); currentnode != dates.back().end(); currentnode++) {
				double M = currentnode->x*e;			
				int k = boost::math::iround(M/delta_x);
				double R = (M-k*delta_x)/V; // eta/V
				for (int l=-1; l<=+1; l++) {
					Transition newtransition((1+R*R)/6.0+((0==l)?(1-R*R):l*R/sqrt(3.0))/2.0);
					auto existingnode = newdate.begin();
					while (existingnode != newdate.end() && newtransition.destination == NULL) {
						if (existingnode->relative_position == k+l) newtransition.destination = &(*existingnode);
						existingnode++;
					}
					if (newtransition.destination == NULL) {
						newdate.push_back(Node(k+l, (k+l)*delta_x));
						newtransition.destination=&newdate.back();
					}
					currentnode->transitions[1+l]=newtransition;
					newtransition.destination->q+=currentnode->q*newtransition.probability*exp(-currentnode->r()*delta_t);
				}
			}
			dates.push_back(newdate);
			double alpha = 0;
			for (auto currentnode = dates.back().begin(); currentnode != dates.back().end(); currentnode++) alpha+=currentnode->q*exp(-currentnode->x*delta_t);
			alpha = log(alpha/(hw_.ratecurve.zerocoupon(*(time+1))))/delta_t;
			std::cout << fixed_percentage(alpha) << "\n";
			for (auto currentnode = dates.back().begin(); currentnode != dates.back().end(); currentnode++) currentnode->alpha = alpha;
			prev_time = *time;
		}
	}
};

std::ostream& operator<< (std::ostream& flux, Tree& tree) {
	for (auto currentdate = tree.dates.begin(); currentdate != tree.dates.end(); currentdate++) {
		for (auto currentnode = currentdate->begin(); currentnode != currentdate->end(); currentnode++) {
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