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
  double x[604], y[604], nx[604], ny[604];
  int n;
  vector<Vector2d> plots;

  std::stringstream ss;
  XFoil *foil;
  virtual void SetUp() {
    n = loadDatFile("sample/CLARK_Y.dat", x, y);
    
    foil = new XFoil();

    if (!foil->initXFoilGeometry(n, x, y, nx, ny)) {
      std::cout << "Initialization error!" << std::endl;
    }
    if (!foil->initXFoilAnalysis(100000, 0, 0.0, 9.0, 1.0, 1.0, 1, 1, true, ss)) {
      std::cout << "Initialization error!" << std::endl;
    }

    for (int i=1; i<=foil->n; i++) {
      Vector2d plot = {foil->x[i], foil->y[i]};
      plots.push_back(plot);
    }
  }  
};

TEST_F(DatGoogleTest, test_cang) {
  //when
  PairIndex actual = foil->cang(plots);

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
  XFoil::C_f actual = foil->cft(6.0, 1E+5, 0.01);

  //then
  ASSERT_DOUBLE_EQ(-0.00021874525090522493, actual.cf);
  ASSERT_DOUBLE_EQ(-2.2177015917117524e-06, actual.hk);
  ASSERT_DOUBLE_EQ(-9.7729505058686228e-13, actual.rt);
  ASSERT_DOUBLE_EQ(2.1840616807413529e-05, actual.msq);
}

//TODO blmidのリファクタリング
TEST_F(DatGoogleTest, test_blmid_turbulent) {
  //given
  foil->hk1 = foil->hk2 = 6.0;
  foil->rt1 = foil->rt2 = 1E+5;
  foil->m1 = foil->m2 = 0.01;
  
  //when
  foil->blmid(1);

  //then
  ASSERT_DOUBLE_EQ(-6.8333333333333339e-07, foil->cfm);
  ASSERT_DOUBLE_EQ(0, foil->cfm_u1);
  ASSERT_DOUBLE_EQ(0, foil->cfm_t1);
  ASSERT_DOUBLE_EQ(0, foil->cfm_d1);

  ASSERT_DOUBLE_EQ(0, foil->cfm_u2);
  ASSERT_DOUBLE_EQ(0, foil->cfm_t2);
  ASSERT_DOUBLE_EQ(0, foil->cfm_d2);

  ASSERT_DOUBLE_EQ(0, foil->cfm_ms);
  ASSERT_DOUBLE_EQ(0, foil->cfm_re);
}

TEST_F(DatGoogleTest, test_blmid_laminar) {
  //given
  foil->hk1 = foil->hk2 = 6.0;
  foil->rt1 = foil->rt2 = 1E+5;
  foil->m1 = foil->m2 = 0.01;
  
  //when
  foil->blmid(2);

  //then
  ASSERT_DOUBLE_EQ(-6.8333333333333339e-07, foil->cfm);
  ASSERT_DOUBLE_EQ(0, foil->cfm_u1);
  ASSERT_DOUBLE_EQ(0, foil->cfm_t1);
  ASSERT_DOUBLE_EQ(0, foil->cfm_d1);

  ASSERT_DOUBLE_EQ(0, foil->cfm_u2);
  ASSERT_DOUBLE_EQ(0, foil->cfm_t2);
  ASSERT_DOUBLE_EQ(0, foil->cfm_d2);

  ASSERT_DOUBLE_EQ(0, foil->cfm_ms);
  ASSERT_DOUBLE_EQ(0, foil->cfm_re);
}

TEST_F(DatGoogleTest, test_blmid_turbulent_wake) {
  //given
  foil->hk1 = foil->hk2 = 6.0;
  foil->rt1 = foil->rt2 = 1E+5;
  foil->m1 = foil->m2 = 0.01;
  
  //when
  foil->blmid(3);

  //then
  ASSERT_DOUBLE_EQ(0, foil->cfm);
  ASSERT_DOUBLE_EQ(0, foil->cfm_u1);
  ASSERT_DOUBLE_EQ(0, foil->cfm_t1);
  ASSERT_DOUBLE_EQ(0, foil->cfm_d1);

  ASSERT_DOUBLE_EQ(0, foil->cfm_u2);
  ASSERT_DOUBLE_EQ(0, foil->cfm_t2);
  ASSERT_DOUBLE_EQ(0, foil->cfm_d2);

  ASSERT_DOUBLE_EQ(0, foil->cfm_ms);
  ASSERT_DOUBLE_EQ(0, foil->cfm_re);
}

TEST_F(DatGoogleTest, test_inside_true) {
  //when
  bool actual = foil->inside(x, y, n, 0.1, 0.05);

  //then
  ASSERT_TRUE(actual);
}

TEST_F(DatGoogleTest, test_inside_false) {
  //when
  bool actual = foil->inside(x, y, n, 0.1, 0.35);

  //then
  ASSERT_FALSE(actual);
}

int main() {
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}