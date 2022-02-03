#ifndef LIBSPANNER_CONSTANTS_H
#define LIBSPANNER_CONSTANTS_H

#include "types.h"

const number_t EPSILON = 0.000001;
const number_t INF = std::numeric_limits<double>::infinity();
const size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

const number_t PI = M_PI;
const number_t PI_OVER_TWO = PI / 2;
const number_t SIX_PI_OVER_SEVEN = 6 * PI / 7;
const number_t FOUR_PI_OVER_SEVEN = 4 * PI / 7;
const number_t TAN30 = tan(PI / 6);
//const double COS30 = cos(PI / 6);
const double COT30 = 1/TAN30;
const number_t PI_OVER_FIVE = PI / 5;
const number_t FOUR_PI_OVER_FIVE = 4 * PI / 5;

#endif //LIBLIBSPANNER_CONSTANTS_H
