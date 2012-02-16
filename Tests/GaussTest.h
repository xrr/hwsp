#include "gtest/gtest.h"
#include "Gauss.h"

using ::testing::TestWithParam;
using ::testing::Values;

typedef Gauss* GaussFactory();
Gauss* GaussFactory1() {return new AbramowitzStegunGauss();}
Gauss* GaussFactory2() {return new BoostGauss();}

class GaussTest : public TestWithParam<GaussFactory*> {
public:
	virtual void SetUp() {gauss_= (*GetParam())();}
	virtual void TearDown() {delete gauss_;}
	virtual ~GaussTest() { delete gauss_; }
protected:
	Gauss* gauss_;
};

TEST_P(GaussTest, Erf) {
	EXPECT_NEAR(-1.000000,gauss_->erf(-5),.000001);
	EXPECT_NEAR(-0.842701,gauss_->erf(-1),.000001);
	EXPECT_NEAR(+0.000000,gauss_->erf(+0),.000001);
	EXPECT_NEAR(+0.842701,gauss_->erf(+1),.000001);
	EXPECT_NEAR(+0.995322,gauss_->erf(+2),.000001);
	EXPECT_NEAR(+1.000000,gauss_->erf(+5),.000001);
}

/*TEST_P(GaussTest, CDF) {
	EXPECT_EQ(gauss_->cdf(1),10);
	EXPECT_EQ(gauss_->cdf(2),20);
}

TEST_P(GaussTest, CDF2) {
	EXPECT_EQ(gauss_->cdf(1),10);
	EXPECT_EQ(gauss_->cdf(2),20);
}*/

//INSTANTIATE_TEST_CASE_P(Impls12, GaussTest,
//	Values(&GaussFactory2));

INSTANTIATE_TEST_CASE_P(Impls12, GaussTest,
	Values(&GaussFactory1, &GaussFactory2));


/* ***************** */

//class GaussTest : public testing::Test {
//protected:
//	AbramowitzStegun Gauss1;
//	BoostGauss Gauss2;
//
//	void CDFTester(Gauss* pG) {
//		EXPECT_EQ(pG->cdf(1),0.5135415);
//		EXPECT_EQ(pG->cdf(1.5),0.65412);
//	}
//
//	void ErfTester(Gauss* pG) {
//		EXPECT_EQ(pG->cdf(1),0.5135415);
//		EXPECT_EQ(pG->cdf(1.5),0.65412);
//	}
//
//};
//
//TEST_F(GaussTest, CDFTest) {
//	CDFTester(&Gauss1);
//	CDFTester(&Gauss2);
//}
//
//TEST_F(GaussTest, ErfTest) {
//	ErfTester(&Gauss1);
//	ErfTester(&Gauss2);
//}