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
  int imax;
  double dmax;
  XFoil *foil = new XFoil();
  bool actual = foil->cang(x,y,n,imax,dmax);
  printf("%d\n%lf", imax,dmax);
  ASSERT_TRUE(actual);
  ASSERT_EQ(78, imax);
  ASSERT_DOUBLE_EQ(9.9742071999702819, dmax);

}

int main() {
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}