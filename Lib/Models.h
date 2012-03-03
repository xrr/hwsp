#pragma once
#include <vector>

class RateCurve {
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
	HullWhite(double _a, double _sigma) : a(_a), sigma(_sigma) {};
};