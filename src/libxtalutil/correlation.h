
#ifndef CORRELATION_H_
#define CORRELATION_H_

#include <math.h>
#include <vector>

/*----------------------------------------------------------------------------*/

float_t stdDeviation(const std::vector<float_t> &x);

float_t coVariance(const std::vector<float_t> &x, const std::vector<float_t> &y);

float_t corrPearson(const std::vector<float_t> &x, const std::vector<float_t> &y);

float_t corrSpears(const std::vector<float_t> &x, const std::vector<float_t> &y);

//void linearRegression(const std::vector<float_t> &x,
//      const std::vector<float_t> &y, float_t &a, float_t &b);

/*----------------------------------------------------------------------------*/

#endif /* CORRELATION_H_ */
