#pragma once

#include <iostream>
#include <sstream>

#include <math.h>

#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/normal.hpp>

#define M_SQRT1_2 7.0710678118654752440E-1 //2^(-1/2)

std::string percentage (double value, int precision = 3) {
	std::ostringstream out (std::ostringstream::out);
	if (0==value) precision--;
	out.precision(precision);
	out << std::showpoint << std::showpos << std::setw(precision) << std::left;
	out << 100*value << "%";
	return out.str();
};


// todo : stop after a given time ?
template <typename S, typename T, typename U>
S dichotomicRootSearch(U f, T f_target, T precision, S min, S max) {
	S target;
	T f_min = f(min)-f_target;
	T f_max = f(max)-f_target;
	T f_target_temp = abs(precision)+123; // arbitrary positive constant
	while (abs(f_target_temp) > precision) { // seeking for function's zero
		target = (min+max) / 2 ;
		f_target_temp = f(target)-f_target;
		if (f_min*f_target_temp < 0) { 
			max = target;
			f_max = f_target_temp;
		} else  {
			min = target;
			f_max = f_target_temp;
		}
	}
	return target;
};


// erf(z) = 2/sqrt(pi) \quad_0^x exp(-t^2) dt
// cdf(z) = 1/sqrt(2*pi) \quad_-inf^x exp(-t^2) dt
class Gauss {
public:
	static double _erf(double){};
	virtual double erf(double)=0;
	virtual double cdf(double x) {return (1+erf(x*M_SQRT1_2))/2;};
};

// erf(0.01) = 0.0112834772 erf(3.7) = 0.9999998325
// Abramowitz/Stegun: p299, |erf(z)-erf| <= 1.5*10^(-7)
class AbramowitzStegunGauss : public Gauss {
public:
	virtual double erf(double x) {return _erf(x);};
	static double _erf(double x) {
		if (x<0) {
			return -_erf(-x);
		} else {
			double y = 1.0 / ( 1.0 + 0.3275911 * x);   
			return 1 - (((((
				+ 1.061405429  * y
				- 1.453152027) * y
				+ 1.421413741) * y
				- 0.284496736) * y 
				+ 0.254829592) * y) 
				* exp (-x * x);
		}
	};
};

// http://boost.org/doc/libs/1_48_0/boost/math/special_functions/erf.hpp
class BoostGauss : public Gauss {
public:
	virtual double erf(double x) {return _erf(x);};
	static double _erf(double x) {return boost::math::erf(x);}
	virtual double cdf(double x) {boost::math::normal norm;return boost::math::cdf(norm, x);};
};