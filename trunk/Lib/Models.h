#pragma once
#include <vector>

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