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
