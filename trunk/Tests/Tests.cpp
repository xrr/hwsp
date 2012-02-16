#include <iostream> 
#include "gtest/gtest.h"

#include "Payoff.h"

TEST(sample_test_case, sample_test) {
	Payoff myPayoff;
	EXPECT_EQ(1, 2-1);
}

int main(int argc, char** argv) {
	testing::InitGoogleTest(&argc, argv);
	RUN_ALL_TESTS();
	std::getchar();
	return 0;
} 