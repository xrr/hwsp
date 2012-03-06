#pragma once
#include <iostream>
#include <sstream>

std::string fixed_percentage (double value, int precision = 3) {
	std::ostringstream out (std::ostringstream::out);
	if (0==value) precision--;
	out.precision(precision);
	out << std::showpoint << std::showpos << std::setw(precision) << std::left;
	out << 100*value << "%";
	return out.str();
};

//Romain test
//dichotomic root search of r*
//template <typename S, typename T, typename U>
//S dichotomicRootSearch(U* pFunc, S target, T precision, S min, S max) {
////S dichotomicRootSearch(S target, T precision, S min, S max) {
//
//		//T f_min = (*pFunc)(swaption, min)-f_target;
//		//T f_max = (*pFunc)(swaption, max)-f_target;
//
//		//T precision = 0.00001;
//		//T f_target_temp = abs(precision)+123;
//		//while (abs(f_target_temp) > precision) { // seeking for function's zero
//		//	target = (min+max) / 2 ;
//		//	f_target_temp = f(swaption, target)-f_target;
//		//	if (f_min*f_target_temp < 0) { 
//		//		max = target;
//		//		f_max = f_target_temp;
//		//	} else  {
//		//		min = target;
//		//		f_max = f_target_temp;
//		//	}
//		//}
//		return target;
//};