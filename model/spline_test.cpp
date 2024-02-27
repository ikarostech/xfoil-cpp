#include <gtest/gtest.h>
#include <vector>

#include "spline.hpp"

TEST(scalc, test_value) {
    //given
    Eigen::MatrixX2d points(4, 2);
    points <<   0.0, 0.0,
                1.0, 0.0,
                2.0, 2.0,
                3.0, 3.0;

    //when
    Eigen::VectorXd actual = spline::scalc(points, 4, 4);

    //then
    ASSERT_EQ(0, actual[0]);
    ASSERT_EQ(1.0, actual[1]);
    ASSERT_EQ(3.2360679774997898, actual[2]);
    ASSERT_EQ(4.6502815398728847, actual[3]);
}

TEST(seval, test_value) {
    //given
    double ss = 1.03347;

    Eigen::VectorXd x {{0,1,0.99344,0.98165,0.96782,0.95254,0.93656,0.92034,0.90394,0.88745,0.87096,0.85444,0.83791,0.82138,0.80484,0.78829,0.77175,0.7552,0.73864,0.72207,0.70549,0.6889,0.67234,0.65579,0.63925,0.62277,0.6063,0.5898,0.57336,0.55689,0.54041,0.52391,0.5074,0.49093,0.47448,0.45804,0.44164,0.42531,0.409,0.3927,0.37637,0.36001,0.3436,0.32714,0.31066,0.29417,0.27783,0.26173,0.24588,0.23032,0.21499,0.1999,0.18499,0.17016,0.15538,0.14076,0.12636,0.11227,0.0986,0.08537,0.07278,0.06126,0.05105,0.04218,0.0346,0.02824,0.02293,0.01848,0.01473,0.01155,0.00885,0.0066,0.00472,0.00318,0.00194,0.001,0.00035,-0,-4e-05,0.00022,0.00082,0.00173,0.00293,0.00443,0.00623,0.00837,0.01086,0.01379,0.01726,0.02137,0.02632,0.03239,0.03973,0.04842,0.05871,0.07083,0.08476,0.09979,0.11534,0.13129,0.14744,0.16378,0.18028,0.19687,0.21359,0.23047,0.24741,0.26443,0.28157,0.29876,0.31597,0.33319,0.35042,0.36766,0.38492,0.40218,0.41944,0.43669,0.45393,0.47116,0.48838,0.50557,0.52274,0.5399,0.55706,0.57423,0.59142,0.60864,0.62583,0.64299,0.66015,0.67729,0.69444,0.71158,0.72873,0.74587,0.76302,0.78018,0.79735,0.81454,0.83176,0.84894,0.8661,0.88324,0.90035,0.91741,0.93432,0.95093,0.9668,0.98113,0.99325,1}};
    Eigen::VectorXd xs {{-0,-0.972634,-0.973106,-0.97329,-0.973675,-0.974506,-0.97517,-0.976281,-0.97713,-0.977873,-0.978807,-0.979594,-0.980407,-0.981287,-0.982125,-0.98301,-0.983839,-0.98472,-0.985442,-0.986322,-0.987174,-0.987855,-0.988636,-0.989521,-0.990338,-0.991229,-0.992151,-0.992939,-0.993808,-0.994625,-0.995399,-0.996039,-0.996595,-0.997245,-0.99782,-0.998302,-0.998795,-0.999218,-0.999574,-0.999815,-0.999964,-1,-0.999963,-0.999843,-0.999683,-0.999492,-0.999178,-0.998675,-0.997801,-0.996514,-0.994722,-0.992093,-0.98887,-0.985078,-0.980996,-0.976318,-0.970699,-0.964098,-0.955861,-0.947104,-0.937009,-0.923422,-0.904611,-0.881308,-0.852305,-0.817688,-0.779868,-0.745129,-0.720148,-0.707049,-0.687086,-0.656067,-0.611032,-0.552088,-0.477249,-0.379435,-0.259228,-0.105845,0.0586071,0.226305,0.38403,0.499802,0.594212,0.67308,0.740235,0.7961,0.848922,0.901883,0.93924,0.964494,0.977229,0.98165,0.9847,0.990812,0.995725,0.998031,0.998816,0.999263,0.999714,0.999917,1,0.99998,0.999892,0.999752,0.999631,0.999516,0.999403,0.99931,0.999276,0.999287,0.999295,0.999324,0.999335,0.999333,0.999335,0.999334,0.999334,0.999332,0.999336,0.999319,0.999316,0.999337,0.999308,0.99933,0.999332,0.999305,0.999353,0.999329,0.999311,0.999331,0.999325,0.999326,0.999326,0.999326,0.999326,0.999325,0.999326,0.999327,0.999328,0.999331,0.99933,0.999327,0.999327,0.99932,0.999336,0.999324,0.99931,0.999337,0.999319,0.999332,0.999331,0.999298}};
    Eigen::VectorXd s {{-0,0,0.00674294,0.0188563,0.0330653,0.04875,0.0651438,0.0817672,0.0985577,0.115428,0.132283,0.149153,0.166021,0.182873,0.199722,0.216565,0.233384,0.250199,0.267009,0.283817,0.300619,0.317419,0.334176,0.350909,0.367617,0.384251,0.400858,0.417482,0.434032,0.450598,0.46716,0.483731,0.500302,0.516823,0.533313,0.549786,0.566209,0.582555,0.598875,0.61518,0.631512,0.647872,0.664282,0.680743,0.697227,0.713724,0.730075,0.746192,0.762069,0.777673,0.793069,0.808258,0.823311,0.838336,0.853371,0.868308,0.883099,0.897662,0.9119,0.925805,0.939165,0.951543,0.962705,0.972633,0.981368,0.988979,0.995621,1.00146,1.00658,1.01104,1.01491,1.01825,1.02121,1.02385,1.02625,1.02843,1.03045,1.03232,1.03414,1.03599,1.03793,1.03997,1.04216,1.04452,1.04706,1.04984,1.05287,1.05621,1.05998,1.06429,1.06938,1.07557,1.08304,1.09184,1.10219,1.11435,1.1283,1.14335,1.1589,1.17486,1.19101,1.20735,1.22385,1.24044,1.25717,1.27405,1.291,1.30803,1.32519,1.34239,1.35961,1.37684,1.39408,1.41133,1.42861,1.44588,1.46315,1.48041,1.49766,1.5149,1.53214,1.54934,1.56652,1.58369,1.60086,1.61804,1.63525,1.65248,1.66968,1.68685,1.70402,1.72117,1.73833,1.75549,1.77265,1.7898,1.80696,1.82413,1.84131,1.85852,1.87575,1.89294,1.91011,1.92726,1.94438,1.96145,1.97838,1.995,2.01088,2.02522,2.03735,2.0441}};
        
    //when
    double actual = spline::seval(ss, x, xs, s, x.size());
    //then
    ASSERT_EQ(-5.9902044244614357e-05, actual);
}

TEST(deval, test_value) {
    //given
    double ss = 1.03347;

    Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(new double[152]{0,1,0.99344,0.98165,0.96782,0.95254,0.93656,0.92034,0.90394,0.88745,0.87096,0.85444,0.83791,0.82138,0.80484,0.78829,0.77175,0.7552,0.73864,0.72207,0.70549,0.6889,0.67234,0.65579,0.63925,0.62277,0.6063,0.5898,0.57336,0.55689,0.54041,0.52391,0.5074,0.49093,0.47448,0.45804,0.44164,0.42531,0.409,0.3927,0.37637,0.36001,0.3436,0.32714,0.31066,0.29417,0.27783,0.26173,0.24588,0.23032,0.21499,0.1999,0.18499,0.17016,0.15538,0.14076,0.12636,0.11227,0.0986,0.08537,0.07278,0.06126,0.05105,0.04218,0.0346,0.02824,0.02293,0.01848,0.01473,0.01155,0.00885,0.0066,0.00472,0.00318,0.00194,0.001,0.00035,-0,-4e-05,0.00022,0.00082,0.00173,0.00293,0.00443,0.00623,0.00837,0.01086,0.01379,0.01726,0.02137,0.02632,0.03239,0.03973,0.04842,0.05871,0.07083,0.08476,0.09979,0.11534,0.13129,0.14744,0.16378,0.18028,0.19687,0.21359,0.23047,0.24741,0.26443,0.28157,0.29876,0.31597,0.33319,0.35042,0.36766,0.38492,0.40218,0.41944,0.43669,0.45393,0.47116,0.48838,0.50557,0.52274,0.5399,0.55706,0.57423,0.59142,0.60864,0.62583,0.64299,0.66015,0.67729,0.69444,0.71158,0.72873,0.74587,0.76302,0.78018,0.79735,0.81454,0.83176,0.84894,0.8661,0.88324,0.90035,0.91741,0.93432,0.95093,0.9668,0.98113,0.99325,1} , 152);
    Eigen::VectorXd xs = Eigen::Map<Eigen::VectorXd>(new double[152]{-0,-0.972634,-0.973106,-0.97329,-0.973675,-0.974506,-0.97517,-0.976281,-0.97713,-0.977873,-0.978807,-0.979594,-0.980407,-0.981287,-0.982125,-0.98301,-0.983839,-0.98472,-0.985442,-0.986322,-0.987174,-0.987855,-0.988636,-0.989521,-0.990338,-0.991229,-0.992151,-0.992939,-0.993808,-0.994625,-0.995399,-0.996039,-0.996595,-0.997245,-0.99782,-0.998302,-0.998795,-0.999218,-0.999574,-0.999815,-0.999964,-1,-0.999963,-0.999843,-0.999683,-0.999492,-0.999178,-0.998675,-0.997801,-0.996514,-0.994722,-0.992093,-0.98887,-0.985078,-0.980996,-0.976318,-0.970699,-0.964098,-0.955861,-0.947104,-0.937009,-0.923422,-0.904611,-0.881308,-0.852305,-0.817688,-0.779868,-0.745129,-0.720148,-0.707049,-0.687086,-0.656067,-0.611032,-0.552088,-0.477249,-0.379435,-0.259228,-0.105845,0.0586071,0.226305,0.38403,0.499802,0.594212,0.67308,0.740235,0.7961,0.848922,0.901883,0.93924,0.964494,0.977229,0.98165,0.9847,0.990812,0.995725,0.998031,0.998816,0.999263,0.999714,0.999917,1,0.99998,0.999892,0.999752,0.999631,0.999516,0.999403,0.99931,0.999276,0.999287,0.999295,0.999324,0.999335,0.999333,0.999335,0.999334,0.999334,0.999332,0.999336,0.999319,0.999316,0.999337,0.999308,0.99933,0.999332,0.999305,0.999353,0.999329,0.999311,0.999331,0.999325,0.999326,0.999326,0.999326,0.999326,0.999325,0.999326,0.999327,0.999328,0.999331,0.99933,0.999327,0.999327,0.99932,0.999336,0.999324,0.99931,0.999337,0.999319,0.999332,0.999331,0.999298}, 152);
    Eigen::VectorXd s = Eigen::Map<Eigen::VectorXd>(new double[152]{-0,0,0.00674294,0.0188563,0.0330653,0.04875,0.0651438,0.0817672,0.0985577,0.115428,0.132283,0.149153,0.166021,0.182873,0.199722,0.216565,0.233384,0.250199,0.267009,0.283817,0.300619,0.317419,0.334176,0.350909,0.367617,0.384251,0.400858,0.417482,0.434032,0.450598,0.46716,0.483731,0.500302,0.516823,0.533313,0.549786,0.566209,0.582555,0.598875,0.61518,0.631512,0.647872,0.664282,0.680743,0.697227,0.713724,0.730075,0.746192,0.762069,0.777673,0.793069,0.808258,0.823311,0.838336,0.853371,0.868308,0.883099,0.897662,0.9119,0.925805,0.939165,0.951543,0.962705,0.972633,0.981368,0.988979,0.995621,1.00146,1.00658,1.01104,1.01491,1.01825,1.02121,1.02385,1.02625,1.02843,1.03045,1.03232,1.03414,1.03599,1.03793,1.03997,1.04216,1.04452,1.04706,1.04984,1.05287,1.05621,1.05998,1.06429,1.06938,1.07557,1.08304,1.09184,1.10219,1.11435,1.1283,1.14335,1.1589,1.17486,1.19101,1.20735,1.22385,1.24044,1.25717,1.27405,1.291,1.30803,1.32519,1.34239,1.35961,1.37684,1.39408,1.41133,1.42861,1.44588,1.46315,1.48041,1.49766,1.5149,1.53214,1.54934,1.56652,1.58369,1.60086,1.61804,1.63525,1.65248,1.66968,1.68685,1.70402,1.72117,1.73833,1.75549,1.77265,1.7898,1.80696,1.82413,1.84131,1.85852,1.87575,1.89294,1.91011,1.92726,1.94438,1.96145,1.97838,1.995,2.01088,2.02522,2.03735,2.0441}, 152);
        
    //when
    double actual = spline::deval(ss, x, xs, s, x.size());
    //then
    ASSERT_EQ(0.00035722664703689773, actual);
}

TEST(d2val, test_value) {
    //given
    double ss = 1.03347;

    Eigen::VectorXd x {{0,1,0.99344,0.98165,0.96782,0.95254,0.93656,0.92034,0.90394,0.88745,0.87096,0.85444,0.83791,0.82138,0.80484,0.78829,0.77175,0.7552,0.73864,0.72207,0.70549,0.6889,0.67234,0.65579,0.63925,0.62277,0.6063,0.5898,0.57336,0.55689,0.54041,0.52391,0.5074,0.49093,0.47448,0.45804,0.44164,0.42531,0.409,0.3927,0.37637,0.36001,0.3436,0.32714,0.31066,0.29417,0.27783,0.26173,0.24588,0.23032,0.21499,0.1999,0.18499,0.17016,0.15538,0.14076,0.12636,0.11227,0.0986,0.08537,0.07278,0.06126,0.05105,0.04218,0.0346,0.02824,0.02293,0.01848,0.01473,0.01155,0.00885,0.0066,0.00472,0.00318,0.00194,0.001,0.00035,-0,-4e-05,0.00022,0.00082,0.00173,0.00293,0.00443,0.00623,0.00837,0.01086,0.01379,0.01726,0.02137,0.02632,0.03239,0.03973,0.04842,0.05871,0.07083,0.08476,0.09979,0.11534,0.13129,0.14744,0.16378,0.18028,0.19687,0.21359,0.23047,0.24741,0.26443,0.28157,0.29876,0.31597,0.33319,0.35042,0.36766,0.38492,0.40218,0.41944,0.43669,0.45393,0.47116,0.48838,0.50557,0.52274,0.5399,0.55706,0.57423,0.59142,0.60864,0.62583,0.64299,0.66015,0.67729,0.69444,0.71158,0.72873,0.74587,0.76302,0.78018,0.79735,0.81454,0.83176,0.84894,0.8661,0.88324,0.90035,0.91741,0.93432,0.95093,0.9668,0.98113,0.99325,1}};
    Eigen::VectorXd xs {{-0,-0.972634,-0.973106,-0.97329,-0.973675,-0.974506,-0.97517,-0.976281,-0.97713,-0.977873,-0.978807,-0.979594,-0.980407,-0.981287,-0.982125,-0.98301,-0.983839,-0.98472,-0.985442,-0.986322,-0.987174,-0.987855,-0.988636,-0.989521,-0.990338,-0.991229,-0.992151,-0.992939,-0.993808,-0.994625,-0.995399,-0.996039,-0.996595,-0.997245,-0.99782,-0.998302,-0.998795,-0.999218,-0.999574,-0.999815,-0.999964,-1,-0.999963,-0.999843,-0.999683,-0.999492,-0.999178,-0.998675,-0.997801,-0.996514,-0.994722,-0.992093,-0.98887,-0.985078,-0.980996,-0.976318,-0.970699,-0.964098,-0.955861,-0.947104,-0.937009,-0.923422,-0.904611,-0.881308,-0.852305,-0.817688,-0.779868,-0.745129,-0.720148,-0.707049,-0.687086,-0.656067,-0.611032,-0.552088,-0.477249,-0.379435,-0.259228,-0.105845,0.0586071,0.226305,0.38403,0.499802,0.594212,0.67308,0.740235,0.7961,0.848922,0.901883,0.93924,0.964494,0.977229,0.98165,0.9847,0.990812,0.995725,0.998031,0.998816,0.999263,0.999714,0.999917,1,0.99998,0.999892,0.999752,0.999631,0.999516,0.999403,0.99931,0.999276,0.999287,0.999295,0.999324,0.999335,0.999333,0.999335,0.999334,0.999334,0.999332,0.999336,0.999319,0.999316,0.999337,0.999308,0.99933,0.999332,0.999305,0.999353,0.999329,0.999311,0.999331,0.999325,0.999326,0.999326,0.999326,0.999326,0.999325,0.999326,0.999327,0.999328,0.999331,0.99933,0.999327,0.999327,0.99932,0.999336,0.999324,0.99931,0.999337,0.999319,0.999332,0.999331,0.999298}};
    Eigen::VectorXd s {{-0,0,0.00674294,0.0188563,0.0330653,0.04875,0.0651438,0.0817672,0.0985577,0.115428,0.132283,0.149153,0.166021,0.182873,0.199722,0.216565,0.233384,0.250199,0.267009,0.283817,0.300619,0.317419,0.334176,0.350909,0.367617,0.384251,0.400858,0.417482,0.434032,0.450598,0.46716,0.483731,0.500302,0.516823,0.533313,0.549786,0.566209,0.582555,0.598875,0.61518,0.631512,0.647872,0.664282,0.680743,0.697227,0.713724,0.730075,0.746192,0.762069,0.777673,0.793069,0.808258,0.823311,0.838336,0.853371,0.868308,0.883099,0.897662,0.9119,0.925805,0.939165,0.951543,0.962705,0.972633,0.981368,0.988979,0.995621,1.00146,1.00658,1.01104,1.01491,1.01825,1.02121,1.02385,1.02625,1.02843,1.03045,1.03232,1.03414,1.03599,1.03793,1.03997,1.04216,1.04452,1.04706,1.04984,1.05287,1.05621,1.05998,1.06429,1.06938,1.07557,1.08304,1.09184,1.10219,1.11435,1.1283,1.14335,1.1589,1.17486,1.19101,1.20735,1.22385,1.24044,1.25717,1.27405,1.291,1.30803,1.32519,1.34239,1.35961,1.37684,1.39408,1.41133,1.42861,1.44588,1.46315,1.48041,1.49766,1.5149,1.53214,1.54934,1.56652,1.58369,1.60086,1.61804,1.63525,1.65248,1.66968,1.68685,1.70402,1.72117,1.73833,1.75549,1.77265,1.7898,1.80696,1.82413,1.84131,1.85852,1.87575,1.89294,1.91011,1.92726,1.94438,1.96145,1.97838,1.995,2.01088,2.02522,2.03735,2.0441}};
        
    //when
    double actual = spline::d2val(ss, x, xs, s, x.size());
    //then
    ASSERT_EQ(88.931575080506789, actual);
}

TEST(sinvrt, test_value) {
    //given
    Eigen::VectorXd x {{0,1,0.99344,0.98165,0.96782,0.95254,0.93656,0.92034,0.90394,0.88745,0.87096,0.85444,0.83791,0.82138,0.80484,0.78829,0.77175,0.7552,0.73864,0.72207,0.70549,0.6889,0.67234,0.65579,0.63925,0.62277,0.6063,0.5898,0.57336,0.55689,0.54041,0.52391,0.5074,0.49093,0.47448,0.45804,0.44164,0.42531,0.409,0.3927,0.37637,0.36001,0.3436,0.32714,0.31066,0.29417,0.27783,0.26173,0.24588,0.23032,0.21499,0.1999,0.18499,0.17016,0.15538,0.14076,0.12636,0.11227,0.0986,0.08537,0.07278,0.06126,0.05105,0.04218,0.0346,0.02824,0.02293,0.01848,0.01473,0.01155,0.00885,0.0066,0.00472,0.00318,0.00194,0.001,0.00035,-0,-4e-05,0.00022,0.00082,0.00173,0.00293,0.00443,0.00623,0.00837,0.01086,0.01379,0.01726,0.02137,0.02632,0.03239,0.03973,0.04842,0.05871,0.07083,0.08476,0.09979,0.11534,0.13129,0.14744,0.16378,0.18028,0.19687,0.21359,0.23047,0.24741,0.26443,0.28157,0.29876,0.31597,0.33319,0.35042,0.36766,0.38492,0.40218,0.41944,0.43669,0.45393,0.47116,0.48838,0.50557,0.52274,0.5399,0.55706,0.57423,0.59142,0.60864,0.62583,0.64299,0.66015,0.67729,0.69444,0.71158,0.72873,0.74587,0.76302,0.78018,0.79735,0.81454,0.83176,0.84894,0.8661,0.88324,0.90035,0.91741,0.93432,0.95093,0.9668,0.98113,0.99325,1}};
    Eigen::VectorXd xs {{-0,-0.972634,-0.973106,-0.97329,-0.973675,-0.974506,-0.97517,-0.976281,-0.97713,-0.977873,-0.978807,-0.979594,-0.980407,-0.981287,-0.982125,-0.98301,-0.983839,-0.98472,-0.985442,-0.986322,-0.987174,-0.987855,-0.988636,-0.989521,-0.990338,-0.991229,-0.992151,-0.992939,-0.993808,-0.994625,-0.995399,-0.996039,-0.996595,-0.997245,-0.99782,-0.998302,-0.998795,-0.999218,-0.999574,-0.999815,-0.999964,-1,-0.999963,-0.999843,-0.999683,-0.999492,-0.999178,-0.998675,-0.997801,-0.996514,-0.994722,-0.992093,-0.98887,-0.985078,-0.980996,-0.976318,-0.970699,-0.964098,-0.955861,-0.947104,-0.937009,-0.923422,-0.904611,-0.881308,-0.852305,-0.817688,-0.779868,-0.745129,-0.720148,-0.707049,-0.687086,-0.656067,-0.611032,-0.552088,-0.477249,-0.379435,-0.259228,-0.105845,0.0586071,0.226305,0.38403,0.499802,0.594212,0.67308,0.740235,0.7961,0.848922,0.901883,0.93924,0.964494,0.977229,0.98165,0.9847,0.990812,0.995725,0.998031,0.998816,0.999263,0.999714,0.999917,1,0.99998,0.999892,0.999752,0.999631,0.999516,0.999403,0.99931,0.999276,0.999287,0.999295,0.999324,0.999335,0.999333,0.999335,0.999334,0.999334,0.999332,0.999336,0.999319,0.999316,0.999337,0.999308,0.99933,0.999332,0.999305,0.999353,0.999329,0.999311,0.999331,0.999325,0.999326,0.999326,0.999326,0.999326,0.999325,0.999326,0.999327,0.999328,0.999331,0.99933,0.999327,0.999327,0.99932,0.999336,0.999324,0.99931,0.999337,0.999319,0.999332,0.999331,0.999298}};
    Eigen::VectorXd s {{-0,0,0.00674294,0.0188563,0.0330653,0.04875,0.0651438,0.0817672,0.0985577,0.115428,0.132283,0.149153,0.166021,0.182873,0.199722,0.216565,0.233384,0.250199,0.267009,0.283817,0.300619,0.317419,0.334176,0.350909,0.367617,0.384251,0.400858,0.417482,0.434032,0.450598,0.46716,0.483731,0.500302,0.516823,0.533313,0.549786,0.566209,0.582555,0.598875,0.61518,0.631512,0.647872,0.664282,0.680743,0.697227,0.713724,0.730075,0.746192,0.762069,0.777673,0.793069,0.808258,0.823311,0.838336,0.853371,0.868308,0.883099,0.897662,0.9119,0.925805,0.939165,0.951543,0.962705,0.972633,0.981368,0.988979,0.995621,1.00146,1.00658,1.01104,1.01491,1.01825,1.02121,1.02385,1.02625,1.02843,1.03045,1.03232,1.03414,1.03599,1.03793,1.03997,1.04216,1.04452,1.04706,1.04984,1.05287,1.05621,1.05998,1.06429,1.06938,1.07557,1.08304,1.09184,1.10219,1.11435,1.1283,1.14335,1.1589,1.17486,1.19101,1.20735,1.22385,1.24044,1.25717,1.27405,1.291,1.30803,1.32519,1.34239,1.35961,1.37684,1.39408,1.41133,1.42861,1.44588,1.46315,1.48041,1.49766,1.5149,1.53214,1.54934,1.56652,1.58369,1.60086,1.61804,1.63525,1.65248,1.66968,1.68685,1.70402,1.72117,1.73833,1.75549,1.77265,1.7898,1.80696,1.82413,1.84131,1.85852,1.87575,1.89294,1.91011,1.92726,1.94438,1.96145,1.97838,1.995,2.01088,2.02522,2.03735,2.0441}};
        

    //when
    double si;
    spline::sinvrt(si, 0.5, x, xs, s, x.size() - 1);

    //then
    ASSERT_EQ(88.931575080506789, si);
}

TEST(splind, test_case) {
    //TODO xsはxより1つ要素が多いのでそれに対応する
    //given
    Eigen::VectorXd x {{0,1,0.99344,0.98165,0.96782,0.95254,0.93656,0.92034,0.90394,0.88745,0.87096,0.85444,0.83791,0.82138,0.80484,0.78829,0.77175,0.7552,0.73864,0.72207,0.70549,0.6889,0.67234,0.65579,0.63925,0.62277,0.6063,0.5898,0.57336,0.55689,0.54041,0.52391,0.5074,0.49093,0.47448,0.45804,0.44164,0.42531,0.409,0.3927,0.37637,0.36001,0.3436,0.32714,0.31066,0.29417,0.27783,0.26173,0.24588,0.23032,0.21499,0.1999,0.18499,0.17016,0.15538,0.14076,0.12636,0.11227,0.0986,0.08537,0.07278,0.06126,0.05105,0.04218,0.0346,0.02824,0.02293,0.01848,0.01473,0.01155,0.00885,0.0066,0.00472,0.00318,0.00194,0.001,0.00035,-0,-4e-05,0.00022,0.00082,0.00173,0.00293,0.00443,0.00623,0.00837,0.01086,0.01379,0.01726,0.02137,0.02632,0.03239,0.03973,0.04842,0.05871,0.07083,0.08476,0.09979,0.11534,0.13129,0.14744,0.16378,0.18028,0.19687,0.21359,0.23047,0.24741,0.26443,0.28157,0.29876,0.31597,0.33319,0.35042,0.36766,0.38492,0.40218,0.41944,0.43669,0.45393,0.47116,0.48838,0.50557,0.52274,0.5399,0.55706,0.57423,0.59142,0.60864,0.62583,0.64299,0.66015,0.67729,0.69444,0.71158,0.72873,0.74587,0.76302,0.78018,0.79735,0.81454,0.83176,0.84894,0.8661,0.88324,0.90035,0.91741,0.93432,0.95093,0.9668,0.98113,0.99325,1}};
    Eigen::VectorXd s {{-0,0,0.00674294,0.0188563,0.0330653,0.04875,0.0651438,0.0817672,0.0985577,0.115428,0.132283,0.149153,0.166021,0.182873,0.199722,0.216565,0.233384,0.250199,0.267009,0.283817,0.300619,0.317419,0.334176,0.350909,0.367617,0.384251,0.400858,0.417482,0.434032,0.450598,0.46716,0.483731,0.500302,0.516823,0.533313,0.549786,0.566209,0.582555,0.598875,0.61518,0.631512,0.647872,0.664282,0.680743,0.697227,0.713724,0.730075,0.746192,0.762069,0.777673,0.793069,0.808258,0.823311,0.838336,0.853371,0.868308,0.883099,0.897662,0.9119,0.925805,0.939165,0.951543,0.962705,0.972633,0.981368,0.988979,0.995621,1.00146,1.00658,1.01104,1.01491,1.01825,1.02121,1.02385,1.02625,1.02843,1.03045,1.03232,1.03414,1.03599,1.03793,1.03997,1.04216,1.04452,1.04706,1.04984,1.05287,1.05621,1.05998,1.06429,1.06938,1.07557,1.08304,1.09184,1.10219,1.11435,1.1283,1.14335,1.1589,1.17486,1.19101,1.20735,1.22385,1.24044,1.25717,1.27405,1.291,1.30803,1.32519,1.34239,1.35961,1.37684,1.39408,1.41133,1.42861,1.44588,1.46315,1.48041,1.49766,1.5149,1.53214,1.54934,1.56652,1.58369,1.60086,1.61804,1.63525,1.65248,1.66968,1.68685,1.70402,1.72117,1.73833,1.75549,1.77265,1.7898,1.80696,1.82413,1.84131,1.85852,1.87575,1.89294,1.91011,1.92726,1.94438,1.96145,1.97838,1.995,2.01088,2.02522,2.03735,2.0441}};
    
    //when
    Eigen::VectorXd xs = spline::splind(x, s, x.size() - spline::INDEX_START_WITH);

    //then
    ASSERT_EQ(0, xs[0]);
    ASSERT_FLOAT_EQ(-0.97263533, xs[1]);
    ASSERT_FLOAT_EQ(-0.97310346, xs[2]);
    ASSERT_FLOAT_EQ(-0.93704528, xs[60]);
    ASSERT_FLOAT_EQ(1.000417, xs[x.size() - 1]);
}

int main() {
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}