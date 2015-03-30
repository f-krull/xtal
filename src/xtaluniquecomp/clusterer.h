#ifndef CLUSTERER_H_
#define CLUSTERER_H_

#include "distmatrix.h"
#include <vector>

/*----------------------------------------------------------------------------*/

template<class T>
class Pairing {
public:
   uint32_t i;
   uint32_t j;
   T distance;
   Pairing(uint32_t i, uint32_t j, T distance);
};

/*----------------------------------------------------------------------------*/

template<class T>
class Clusterer {
private:
protected:
public:
   virtual bool linkNext(uint32_t *i, uint32_t *j, T *dist) = 0;
   virtual ~Clusterer() {}
   static uint32_t clusterAll(DistanceMatrix<T> *dm,
         std::vector<Pairing<T> > *pairs, T maxDist, uint32_t maxSteps);
   static uint32_t clusterAll(DistanceMatrix<T> *dm,
         std::vector<Pairing<T> > *pairs);
};

#endif /*CLUSTERER_H_*/
