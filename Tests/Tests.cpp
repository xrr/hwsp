#include <iostream> 
#include "gtest/gtest.h"

#include "GaussTest.h"

int main(int argc, char** argv) {
	testing::InitGoogleTest(&argc, argv);
	RUN_ALL_TESTS();
	std::getchar();
	return 0;
}