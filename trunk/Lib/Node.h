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

class HullWhite {
public:
	double a, sigma;
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
	double x, r; // displaced, normal rates
	Transition transitions[3]; // down, middle, up
	Node(int _relative_position, double _x) : relative_position(_relative_position), x(_x), r(NULL) {};
};

class Tree {

private:
	HullWhite hw_;
	std::vector<std::vector<Node>> dates;

public:
	std::vector<double> timesteps;

	friend std::ostream& operator<< (std::ostream&, Tree&);

	Tree (std::vector<double> _timesteps, HullWhite hw) : hw_(hw), timesteps(_timesteps) {
		std::vector<Node> date0;
		date0.push_back(Node(0,0));
		dates.push_back(date0);
	};

	void construct(void) {
		for (std::vector<double>::iterator timestep = timesteps.begin(); timestep != timesteps.end(); timestep++) {
			std::vector<Node> newdate;
			double e = exp(-hw_.a*(*timestep));
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
				}
			}
			dates.push_back(newdate);
		}
	}
};

std::ostream& operator<< (std::ostream& flux, Tree& tree) {
	for (auto currentdate = tree.dates.begin(); currentdate != tree.dates.end(); currentdate++) {
		for (auto currentnode = currentdate->begin(); currentnode != currentdate->end(); currentnode++) {
			flux << fixed_percentage(currentnode->x) << " ";
		}
		flux << "\n";
	}
	return flux;
}