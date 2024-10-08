#include <cstdlib>
#include <iostream>
#include "Data3Ddescriptors.h"

using namespace std;

// based on Dragon 6 values provided by kobe Team
double Data3Ddescriptors::mw[110] = {
    0.0841,  0.3333,  0.5779,  0.7504,  0.9001,  1,       1.1665,  1.3322,
    1.582,   1.6803,  1.9143,  2.0237,  2.2465,  2.3389,  2.5787,  2.6703,
    2.9517,  3.3262,  3.2555,  3.3371,  3.7432,  3.9856,  4.2416,  4.3297,
    4.5745,  4.6503,  4.9067,  4.8868,  5.2914,  5.4446,  5.8052,  6.0458,
    6.2381,  6.5745,  6.6528,  6.9775,  7.1164,  7.2956,  7.4027,  7.5957,
    7.7357,  7.9883,  8.1599,  8.4155,  8.5684,  8.8609,  8.9817,  9.3597,
    9.5604,  9.8843,  10.1382, 10.6245, 10.5662, 10.9317, 11.0662, 11.4344,
    11.5659, 11.6666, 11.7326, 12.01,   12.0733, 12.5196, 12.6531, 13.0933,
    13.2327, 13.5304, 13.7327, 13.9267, 14.0661, 14.408,  14.5684, 14.8618,
    15.0664, 15.3072, 15.5043, 15.8393, 16.0047, 16.2431, 16.4005, 16.7019,
    17.0175, 17.2523, 17.4005, 17.4854, 17.4854, 18.4846, 18.5679, 18.8177,
    18.9009, 19.3204, 19.237,  19.8192, 19.7336, 20.3164, 20.2331, 20.5662,
    20.5662, 20.8993, 20.9825, 21.3988, 21.4821, 21.5654, 21.8152, 21.7319,
    21.8152, 22.1482, 21.9817, 22.398,  22.3147, 22.5645};

double Data3Ddescriptors::vdW[110] = {
    0.2634, 0.5583, 1.2269, 0,      0,      1,      0.758,  0.7148, 0.6467,
    0.7434, 2.381,  1.0539, 0,      1.8848, 1.1871, 1.1871, 1.0909, 1.3523,
    4.2328, 0,      0,      0,      0,      0,      0,      0,      0,
    0.8814, 0.5583, 0.5466, 1.3309, 0,      1.2886, 1.396,  1.2886, 1.6778,
    0,      0,      0,      0,      0,      0,      0,      0,      0,
    0.8814, 1.0355, 0.8027, 1.4631, 2.0797, 0,      1.7794, 1.5802, 2.051,
    0,      0,      0,      0,      0,      0,      0,      0,      0,
    0,      0,      0,      0,      0,      0,      0,      0,      0,
    0,      0,      0,      0,      0,      1.0355, 0.931,  0.758,  1.5326,
    1.6778, 0,      0,      0,      0,      0,      0,      0,      0,
    0,      1.3095, 0,      0,      0,      0,      0,      0,      0,
    0,      0,      0,      0,      0,      0,      0,      0,      0,
    0,      0};

double Data3Ddescriptors::neg[110] = {
    0.9418, 0,      0.3236, 0.6582, 0.8291, 1,      1.16,   1.3273, 1.4545,
    1.6364, 0.2036, 0.48,   0.6218, 0.7782, 0.9164, 1.0764, 1.2655, 1.2036,
    0.1636, 0.3455, 0.3709, 0.3964, 0.5055, 0.6036, 0.8,    0.8,    0.9309,
    0.7055, 0.72,   0.8109, 0.88,   0.9527, 1.0255, 1.0945, 1.1709, 1.0582,
    0.1127, 0.2618, 0.2364, 0.3273, 0.5164, 0.4182, 0,      0,      0,
    0,      0.6655, 0.72,   0.7782, 0.8364, 0.8945, 0.9527, 1.0109, 0.8509,
    0.08,   0.2473, 0,      0,      0,      0,      0,      0,      0,
    0,      0,      0,      0,      0,      0,      0,      0,      0,
    0,      0.3564, 0,      0,      0,      0,      0,      0.8,    0.8182,
    0.8327, 0.8509, 0,      0,      0,      0,      0,      0,      0,
    0,      0,      0,      0,      0,      0,      0,      0,      0,
    0,      0,      0,      0,      0,      0,      0,      0,      0,
    0,      0};

double Data3Ddescriptors::pol[110] = {
    0.3807,  0.1136,  13.8068, 3.1818,  1.7216,  1,       0.625,   0.4545,
    0.3182,  0.2216,  13.4091, 6.0227,  3.8636,  3.0568,  2.0625,  1.6477,
    1.2386,  0.9318,  24.6591, 12.9545, 10.1136, 8.2955,  7.0455,  6.5909,
    5.3409,  4.7727,  4.2614,  3.8636,  3.4659,  4.0341,  4.6136,  3.4489,
    2.4489,  2.142,   1.733,   1.4091,  26.875,  15.6818, 12.8977, 10.1705,
    8.9205,  7.2727,  6.4773,  5.4545,  4.8864,  2.7273,  4.0909,  4.0909,
    5.7955,  4.375,   3.75,    3.125,   3.0398,  2.2955,  33.8636, 22.5568,
    17.6705, 16.8182, 16.0227, 17.8409, 17.1023, 16.3636, 15.7386, 13.3523,
    14.4886, 13.9205, 13.4091, 12.8977, 12.3864, 11.9318, 12.4432, 9.2045,
    7.4432,  6.3068,  5.5114,  4.8295,  4.3182,  3.6932,  3.2955,  3.2386,
    4.3182,  3.8636,  4.2045,  3.8636,  3.4091,  3.0114,  27.6705, 21.7614,
    18.2386, 18.2386, 14.4318, 15.5682, 14.0909, 13.9205, 13.2386, 13.0682,
    12.8977, 11.6477, 11.1932, 13.5227, 10.3409, 9.9432,  0,       0,
    0,       0,       0,       0,       0,       0};

double Data3Ddescriptors::ionpol[110] = {
    1.2076, 2.1835, 0.4788, 0.8279, 0.7369, 1.0,    1.2907, 1.2094, 1.5473,
    1.9151, 0.4564, 0.679,  0.5316, 0.7239, 0.9313, 0.92,   1.1516, 1.3996,
    0.3855, 0.5429, 0.5827, 0.6064, 0.5991, 0.6009, 0.6602, 0.7018, 0.6999,
    0.6785, 0.6862, 0.8343, 0.5328, 0.7016, 0.8717, 0.8661, 1.0492, 1.2433,
    0.371,  0.5058, 0.5521, 0.5891, 0.6002, 0.6299, 0.6465, 0.6537, 0.6624,
    0.7404, 0.6728, 0.7987, 0.5139, 0.6522, 0.7645, 0.8001, 0.9282, 1.0772,
    0.3458, 0.4628, 0.4953, 0.4919, 0.486,  0.4907, 0.4929, 0.5012, 0.5036,
    0.5461, 0.5208, 0.5274, 0.5348, 0.5424, 0.5492, 0.5554, 0.4819, 0.6061,
    0.6705, 0.6984, 0.6957, 0.7494, 0.7963, 0.7956, 0.8193, 0.9269, 0.5425,
    0.6587, 0.647,  0.7472, 0,      0.9545, 0.3617, 0.4688, 0.4591, 0.5601,
    0.5231, 0.5501, 0.5564, 0.5352, 0.5305, 0.5321, 0.5504, 0.5579, 0.5701,
    0.5772, 0.5844, 0.5906, 0.4352, 0.5328, 0,      0,      0,      0,
    0,      0};

double Data3Ddescriptors::rcov[110] = {
    0.4079, 0.3684, 1.6842, 1.2632, 1.1053, 1,      0.9342, 0.8684, 0.75,
    0.7632, 2.1842, 1.8553, 1.5921, 1.4605, 1.4079, 1.3816, 1.3421, 1.3947,
    2.6711, 2.3158, 2.2368, 2.1053, 2.0132, 1.8289, 1.8289, 1.7368, 1.6579,
    1.6316, 1.7368, 1.6053, 1.6053, 1.5789, 1.5658, 1.5789, 1.5789, 1.5263,
    2.8947, 2.5658, 2.5,    2.3026, 2.1579, 2.0263, 1.9342, 1.9211, 1.8684,
    1.8289, 1.9079, 1.8947, 1.8684, 1.8289, 1.8289, 1.8158, 1.8289, 1.8421,
    3.2105, 2.8289, 2.7237, 2.6842, 2.6711, 2.6447, 2.6184, 2.6053, 2.6053,
    2.5789, 2.5526, 2.5263, 2.5263, 2.4868, 2.5,    2.4605, 2.4605, 2.3026,
    2.2368, 2.1316, 1.9868, 1.8947, 1.8553, 1.7895, 1.7895, 1.7368, 1.9079,
    1.9211, 1.9474, 1.8421, 1.9737, 1.9737, 3.4211, 2.9079, 2.8289, 2.7105,
    2.6316, 2.5789, 2.5,    2.4605, 2.3684, 2.2237, 0,      0,      0,
    0,      0,      0,      0,      0,      0,      0,      0,      0,
    0,      0};

Data3Ddescriptors::Data3Ddescriptors() {}

double *Data3Ddescriptors::getMW() { return mw; }

double *Data3Ddescriptors::getVDW() { return vdW; }

double *Data3Ddescriptors::getNEG() { return neg; }

double *Data3Ddescriptors::getPOL() { return pol; }

double *Data3Ddescriptors::getIonPOL() { return ionpol; }

double *Data3Ddescriptors::getRCOV() { return rcov; }
