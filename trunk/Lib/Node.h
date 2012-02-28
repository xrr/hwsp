#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <boost/math/special_functions/round.hpp>

class HullWhite {
public:
	double a;
	double sigma;
	//std::vector<double> theta
};

class Node; //double inclusion issue ...

class Transition {
public:
	double probability;
	Node* destination;
};

class Node {
public:
	int relative_position;
	double x; //displaced rate
	double alpha; //displacement
	double r; //rate
	Transition up;
	Transition middle;
	Transition down;
	//todo : initialize alpha,r with NULL or 0 or ??? value
};

class Tree {
private:
	HullWhite hw_;
	std::vector<std::vector<Node>> slices;

public:
	std::vector<double> timesteps;

 Tree (double _timesteps, HullWhite hw) : hw_(hw), timesteps(_timesteps) {
		
		Node node0;
		node0.relative_position=0;
		node0.x=0;
		std::vector<Node> date0;
		date0.push_back(node0);
		slices.push_back(date0);
	};
void construct_1(void) {
		for (std::vector<double>::iterator timestep = timesteps.begin(); timestep != timesteps.end(); timestep++) {
			std::vector<Node> newdate;
			double e = exp(-hw_.a*(*timestep));
			double V = hw_.sigma * sqrt((1-e*e)/(2*hw_.a));
			double delta_x = V * sqrt(3.0);
			for (std::vector<Node>::iterator currentnode = slices.back().begin(); currentnode != slices.back().end(); currentnode++){
				double M = currentnode->x*e;			
				int k = boost::math::iround(M/delta_x);
				for(int l=-1;l<=+1;l++) {
					int relative_position = k+l;
					Transition trans;
					trans.destination = NULL;
					for (std::vector<Node>::iterator existingnode = newdate.begin(); existingnode != newdate.end(); existingnode++) {
						if (existingnode->relative_position == relative_position) {
							trans.destination = &(*existingnode);
						}
					}
					if (trans.destination == NULL) {
						Node creatednode;
						creatednode.relative_position = relative_position;
						creatednode.x=relative_position*delta_x;
						newdate.push_back(creatednode);
						trans.destination=&newdate.back();
					}
					double eta = M-k*delta_x;
					switch (l) {
					case -1:
						trans.probability = 1.0/6 + (eta*eta) / (6*V*V) - eta/(2*sqrt(3.0)*V);
						currentnode->down=trans;
						break;
					case 0:
						trans.probability = 2.0/3 - (eta*eta) / (3*V*V);
						currentnode->middle=trans;
						break;
					case 1:
						trans.probability = 1.0/6 + (eta*eta) / (6*V*V) + eta/(2*sqrt(3.0)*V);
						currentnode->up=trans;
						break;
					} //switch (l)
				} // down, middle, up
			} // currentnode 'j'
			slices.push_back(newdate);
		} // slices.back() (date 'i')
	}; // void construct_1()

	void print(void) {
for (std::vector<std::vector<Node>>::iterator currentslice = slices.begin(); currentslice != slices.end(); currentslice++) {
	for (std::vector<Node>::iterator currentnode = currentslice->begin(); currentnode != currentslice->end(); currentnode++) {
		std::count <<currentnode->x << " ";
	}
	std::count <<std::endl;
}
	}; // void print()


}; //class Tree



