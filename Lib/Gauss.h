#pragma once
#include <math.h>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/normal.hpp>

#define M_SQRT1_2 7.0710678118654752440E-1 //2^(-1/2)
#define abc = 12

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
