#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <boost/math/special_functions/round.hpp>

class HullWhite {
public:
	double a, sigma;
	// std::vector<double> theta
};

class Node; // double inclusion issue ...

class Transition {
public:
	double probability;
	Node *destination;
};

class Node {
public:
	int relative_position;
	double x, r; // displaced, normal rates
	Transition transitions[3]; // down, middle, up
	Node(int _relative_position, double _x) : relative_position(_relative_position), x(_x), r(NULL) {}; // constructor
};

class Tree {
private:
	HullWhite hw_;
	std::vector<std::vector<Node>> slices;

public:
	std::vector<double> timesteps;

	Tree (std::vector<double> _timesteps, HullWhite hw) : hw_(hw), timesteps(_timesteps) {
		std::vector<Node> date0;
		date0.push_back(Node(0,0));
		slices.push_back(date0);
	};

	void construct_1(void) {
		for (std::vector<double>::iterator timestep = timesteps.begin(); timestep != timesteps.end(); timestep++) {
			std::vector<Node> newdate;
			double e = exp(-hw_.a*(*timestep));
			double V = hw_.sigma * sqrt((1-e*e)/(2*hw_.a));
			double delta_x = V * sqrt(3.0);
			for (auto currentnode = slices.back().begin(); currentnode != slices.back().end(); currentnode++) {
				double M = currentnode->x*e;			
				int k = boost::math::iround(M/delta_x);
				double R = (M-k*delta_x)/V; // eta/V
				for(int l=-1;l<=+1;l++) {
					Transition newtransition;
					newtransition.probability = (1+R*R)/6.0+((0==l)?(1-R*R):(0>l?-1:1)*R/sqrt(3.0))/2.0;
					newtransition.destination = NULL; // todo : while loop
					for (auto existingnode = newdate.begin(); existingnode != newdate.end(); existingnode++)
						if (existingnode->relative_position == k+l) newtransition.destination = &(*existingnode);
					if (newtransition.destination == NULL) {
						newdate.push_back(Node(k+l, (k+l)*delta_x));
						newtransition.destination=&newdate.back();
					}
					currentnode->transitions[1+l]=newtransition;
				} // transition l	
			} // currentnode 'j'
			slices.push_back(newdate);
		} //  date 'i' = slices.back()
	}; // void construct_1()
	//todo = check transition pointers ?

	void print(void) {
		std::cout.precision(4);
		//for (auto currentslice = slices.begin(); currentslice != slices.end(); currentslice++) {
		for (auto currentslice = slices.begin(); currentslice != slices.end(); currentslice++) {
			for (auto currentnode = currentslice->begin(); currentnode != currentslice->end(); currentnode++) {
				std::cout << std::showpoint << std::showpos << std::setw(8) << std::left << currentnode->x << " ";
			}
			std::cout << "\n";
		}
	}; // void print()


}; //class Tree



