#pragma once
#include <iostream>
#include <sstream>

std::string percentage (double value, int precision = 3) {
	std::ostringstream out (std::ostringstream::out);
	if (0==value) precision--;
	out.precision(precision);
	out << std::showpoint << std::showpos << std::setw(precision) << std::left;
	out << 100*value << "%";
	return out.str();
};

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