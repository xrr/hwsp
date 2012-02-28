#pragma once
#include <vector>

class Transition {
	public:
		double probability;
		Node* destination;
}

class Node {
public:
	int relative_position;
	double x; //displaced rate
	double alpha; //displacement
	double r; //rate
	Transition up;
	Transition middle;
	Transition down;
};

class Tree {
public:
	Tree(std::vector<double> _timesteps, double HW_a, double HW_sigma,)

	std::vector<double> timesteps;
	

private:
	std::vector<std::vector<Node>>
}

class HullWhite {
public:
	double a;
	double sigma;
	std::vector<double> theta;

}