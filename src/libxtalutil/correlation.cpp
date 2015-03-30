#include "correlation.h"
#include <math.h>
#include <inttypes.h>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <numeric>

/*----------------------------------------------------------------------------*/

static bool cmpSmaller(const std::pair<float_t, uint32_t>& a, const std::pair<
      float_t, uint32_t>& b) {
   return a.first < b.first;
}

/*----------------------------------------------------------------------------*/
/* small values get small ranks. (0.5, 0.2, 0.4) -> (3, 1, 2) */
static std::vector<uint32_t> rank(const std::vector<float_t> &c) {
   std::vector<uint32_t> rank;
   rank.reserve(c.size());

   std::vector<std::pair<float_t, uint32_t> > tmp;
   tmp.resize(c.size());
   for (uint32_t i = 0; i < c.size(); ++i) {
      tmp[i] = (std::pair<float_t, uint32_t>(c[i], i));
   }
   std::sort(tmp.begin(), tmp.end(), cmpSmaller);

   tmp.resize(c.size());
   for (uint32_t i = 0; i < c.size(); ++i) {
      tmp[i] = (std::pair<float_t, uint32_t>(tmp[i].second, i));
   }
   std::sort(tmp.begin(), tmp.end(), cmpSmaller);

   for (uint32_t i = 0; i < c.size(); ++i) {
      rank[i] = tmp[i].second + 1;
   }
   return rank;
}

/*----------------------------------------------------------------------------*/
/* Spearman's rank correlation coefficient */
float_t corrSpears(const std::vector<float_t> &x, const std::vector<float_t> &y) {
   float_t sum = 0;

   assert(x.size() == y.size());
   std::vector<uint32_t> rx = rank(x);
   std::vector<uint32_t> ry = rank(y);

   for (uint32_t i = 0; i < x.size(); ++i) {
      sum += pow(rx[i] - (float_t) ry[i], 2);
   }
   return 1 - ((6 * sum) / (x.size() * (pow(x.size(), 2) - 1)));
}

/*----------------------------------------------------------------------------*/

static float_t avg(const std::vector<float_t> &x) {
   return std::accumulate(x.begin(), x.end(), 0.0f) / x.size();
}

/*----------------------------------------------------------------------------*/

static float_t stdDeviation(const std::vector<float_t> &x, float_t avgX) {
   float_t stdDev = 0;
   //float_t avgX = variance(x);

   for (uint32_t i = 0; i < x.size(); ++i) {
      stdDev += pow(x[i] - avgX, 2);
   }
   return sqrt((stdDev / x.size()));
}

/*----------------------------------------------------------------------------*/

float_t stdDeviation(const std::vector<float_t> &x) {
   return stdDeviation(x, avg(x));
}

/*----------------------------------------------------------------------------*/

static float_t coVariance(const std::vector<float_t> &x, const std::vector<
      float_t> &y, float_t avgX, float_t avgY) {
   float_t cov = 0;

   assert(x.size() == y.size());
   for (uint32_t i = 0; i < x.size(); ++i) {
      cov += (x[i] - avgX) * (y[i] - avgY);
   }
   return cov / x.size();
}

/*----------------------------------------------------------------------------*/

float_t coVariance(const std::vector<float_t> &x, const std::vector<float_t> &y) {
   return coVariance(x, y, avg(x), avg(y));
}

/*----------------------------------------------------------------------------*/
/* Pearson correlation coefficient */
float_t corrPearson(const std::vector<float_t> &x,
      const std::vector<float_t> &y) {
   float_t avgX = avg(x);
   float_t avgY = avg(y);

   float_t numerator = coVariance(x, y, avgX, avgY) ;
   float_t denominator = stdDeviation(x, avgX) * stdDeviation(y, avgY);

   /* dont't divide by zero. return 0 in those cases */
   return denominator != 0.0f ? numerator / denominator : 0;
}


//void linearRegression(const std::vector<float_t> &x,
//      const std::vector<float_t> &y, float_t &a, float_t &b) {
//   float_t avgd = 0;
//   float_t avgn = 0;
//
//   for (uint32_t i = 0; i < x.size(); ++i) {
//      avgd += x[i];
//      avgn += y[i];
//   }
//   avgd /= x.size();
//   avgn /= x.size();
//
//   float_t num = 0;
//   float_t den = 0;
//
//   for (uint32_t i = 0; i < x.size(); ++i) {
//      num += (x[i] - avgd) * (y[i] - avgn);
//      den += pow(x[i] - avgd, 2);
//   }
//
//   a = num / den;
//   b = avgn - a * avgd;
//}
