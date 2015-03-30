#include "clusterer.h"

/*----------------------------------------------------------------------------*/

template <class T>
Pairing<T>::Pairing(uint32_t i, uint32_t j, T distance) {
   this->i = i;
   this->j = j;
   this->distance = distance;
}
