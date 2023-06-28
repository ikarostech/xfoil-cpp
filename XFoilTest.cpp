#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "XFoil.h"

class DatGoogleTest : public ::testing::Test {
protected:
  // TODO テスト用翼型の読み込みは別途Util化する
  virtual int loadDatFile(std::string filename, double x[604], double y[604]) {
    std::ifstream fs(filename);
    if (!fs) {
      std::cout << "Failed to open dat file" << std::endl;
      return -1;
    }

    std::string line;
    std::getline(fs, line);
    std::cout << "Foil name : " << line << std::endl;
    int cnt = 0;
    while (!fs.eof()) {
      std::getline(fs, line);

      line.erase(0, line.find_first_not_of(" \t"));
      int endOfX = line.find_first_of(" \t");
      if (endOfX == -1) continue;

      std::string sx = line.substr(0, endOfX);
      std::string sy = line.substr(endOfX);

      x[cnt] = atof(sx.c_str());
      y[cnt] = atof(sy.c_str());
      cnt++;
    }
    return cnt;
  }
  double x[604], y[604];
  int n;
  virtual void SetUp() {
    n = loadDatFile("sample/CLARK_Y.dat", x, y);
  }  
};

TEST_F(DatGoogleTest, testCang) {
  //given
  XFoil *foil = new XFoil();
  
  //when
  PairIndex actual = foil->cang(x,y,n);

  //then
  ASSERT_EQ(78, actual.index);
  ASSERT_DOUBLE_EQ(9.9742071999702819, actual.value);
}

TEST(cfl_test, hk_over_5_5) {
  //given
  XFoil *foil = new XFoil();

  //when
  XFoil::C_f actual = foil->cfl(6.0, 1E+5);

  //then
  ASSERT_DOUBLE_EQ(-6.8333333333333339e-07, actual.cf);
  ASSERT_DOUBLE_EQ(4.4444444444444448e-08, actual.hk);
  ASSERT_DOUBLE_EQ(6.833333333333334e-12, actual.rt);
  ASSERT_DOUBLE_EQ(0.0, actual.msq);
}

TEST(cfl_test, hk_under_5_5) {
  //given
  XFoil *foil = new XFoil();
  //when
  XFoil::C_f actual = foil->cfl(4.0, 1E+5);

  //then
  ASSERT_DOUBLE_EQ(-2.0927500000000003e-07, actual.cf);
  ASSERT_DOUBLE_EQ(-1.0795950000000001e-06, actual.hk);
  ASSERT_DOUBLE_EQ(2.0927500000000003e-12, actual.rt);
  ASSERT_DOUBLE_EQ(0.0, actual.msq);
}

TEST(cft_test, cal_cft) {
  //given
  XFoil *foil = new XFoil();
  double cf;
  double cf_hk;
  double cf_rt;
  double cf_msq;
                
  //when
  //XFoil::C_f actual = foil->cfl(4.0, 1E+5);
  foil->cft(6.0, 1E+5, 0.01, cf, cf_hk, cf_rt, cf_msq);

  //then
  ASSERT_DOUBLE_EQ(-0.00021874525090522493, cf);
  ASSERT_DOUBLE_EQ(-2.2177015917117524e-06, cf_hk);
  ASSERT_DOUBLE_EQ(-9.7729505058686228e-13, cf_rt);
  ASSERT_DOUBLE_EQ(2.1840616807413529e-05, cf_msq);
}
int main() {
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}