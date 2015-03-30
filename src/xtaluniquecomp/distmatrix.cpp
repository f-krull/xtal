
#ifndef DISTMATRIX_CPP_
#define DISTMATRIX_CPP_

#include "distmatrix.h"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "../libxtalutil/log.h"

/*----------------------------------------------------------------------------*/

template<typename T>
DistanceMatrix<T>::DistanceMatrix(uint32_t numElements, T initValue) {
   this->m_numElements = numElements;
   m_matrix = new T *[numElements - 1];
   for (uint32_t i = 0; i < numElements - 1; i++) {
      m_matrix[i] = new T[i + 1];
      for (uint32_t j = 0; j < i + 1; j++) {
         m_matrix[i][j] = initValue;
      }
   }
   this->m_null = (T) 0;
}

/*----------------------------------------------------------------------------*/

template<typename T>
DistanceMatrix<T>::~DistanceMatrix() {
   for (uint32_t i = 0; i < m_numElements - 1; i++) {
      delete[] m_matrix[i];
   }
   delete[] m_matrix;
}

/*----------------------------------------------------------------------------*/

template<typename T>
T DistanceMatrix<T>::get(uint32_t i, uint32_t j) const {
   assert(i < m_numElements);
   assert(j < m_numElements);
   if (i < j) {
      return m_matrix[j - 1][i];
   } else if (j < i) {
      return m_matrix[i - 1][j];
   } else {
      /* distance from i to i is 0 */
      return (T) 0;
   }
}

/*----------------------------------------------------------------------------*/

template<typename T>
void DistanceMatrix<T>::set(uint32_t i, uint32_t j, T value) {
   assert(i < m_numElements);
   assert(j < m_numElements);
   if (i < j) {
      m_matrix[j - 1][i] = value;
   } else if (j < i) {
      m_matrix[i - 1][j] = value;
   }
}

/*----------------------------------------------------------------------------*/

template<typename T>
void DistanceMatrix<T>::setAll(T value) {
   for (uint32_t i = 0; i < m_numElements; i++) {
      for (uint32_t j = i; j < m_numElements; j++) {
         set(i, j, value);
      }
   }
}

/*----------------------------------------------------------------------------*/

template<typename T>
uint32_t DistanceMatrix<T>::getNumElements() {
   return m_numElements;
}

/*----------------------------------------------------------------------------*/

template<typename T>
void DistanceMatrix<T>::print() {
   for (uint32_t i = 0; i < m_numElements; i++) {
      std::cout << "dm ";
      for (uint32_t j = 0; j < m_numElements; j++) {
         std::cout << std::fixed << std::setprecision(3) << std::setw(7);
         std::cout << get(i, j);
      }
      std::cout << std::endl;
   }
}

/*----------------------------------------------------------------------------*/

template<typename T>
T& DistanceMatrix<T>::operator()(uint32_t i, uint32_t j) {
   assert(i < m_numElements);
   assert(j < m_numElements);
   if (i < j) {
      return m_matrix[j - 1][i];
   } else if (j < i) {
      return m_matrix[i - 1][j];
   } else {
      /* distance from i to i is 0 */
      return m_null;
      /* TODO: i do not like this */
   }
}

/*----------------------------------------------------------------------------*/
template<typename T>
bool DistanceMatrix<T>::save(const char *filename) {
   FILE *f = NULL;

   f = fopen(filename, "wb");
   if (f == NULL) {
      return false;
   }
   bool ok = true;
   ok = ok && fwrite(&m_numElements, sizeof(m_numElements), 1, f) == 1;
   T tmp;
   for (uint32_t i = 0; i < m_numElements; i++) {
      for (uint32_t j = 0; j < m_numElements; j++) {
         tmp = get(i, j);
         ok = ok && fwrite(&tmp, sizeof(T), 1, f) == 1;
      }
   }
   fclose(f);
   return ok;
}

/*----------------------------------------------------------------------------*/
template<typename T>
bool DistanceMatrix<T>::load(const char *filename) {
   FILE *f = NULL;

   f = fopen(filename, "rb");
   if (f == NULL) {
      return false;
   }
   bool ok = true;
   ok = ok && fread(&m_numElements, sizeof(m_numElements), 1, f) == 1;
   T tmp;
   for (uint32_t i = 0; i < m_numElements; i++) {
      for (uint32_t j = 0; j < m_numElements; j++) {
         ok = ok && fread(&tmp, sizeof(T), 1, f) == 1;
         set(i, j, tmp);
      }
   }
   fclose(f);
   return ok;
}

/*----------------------------------------------------------------------------*/
#if 0
template<typename T, typename U>
DistanceMatrix<T>* DistanceMatrixFactory<T, U>::getFilled(
      std::vector<U> *elements, T(*dist)(const U&, const U&)) {
   DistanceMatrix<T>* dm = new DistanceMatrix<T> (elements->size());
   uint32_t numdist;
   uint32_t maxnumdist;

   numdist = 0;
   maxnumdist = (*dm).getNumElements() * (*dm).getNumElements();
   maxnumdist += (*dm).getNumElements();
   maxnumdist /= 2;
   //# pragma omp single
   for (uint32_t i = 0; i < (*dm).getNumElements(); i++) {
      #pragma omp parallel for shared(numdist)
      for (uint32_t j = i + 1; j < (*dm).getNumElements(); j++) {
         (*dm)(i, j) = dist((*elements)[i], (*elements)[j]);

         # pragma omp atomic
         numdist++;

         if (numdist % 1000 == 0) {
            Log::dbg("distance matrix calculated %u/%u distances", numdist,
                  maxnumdist);
         }
      }
   }
   return dm;
}
#endif

template<class T> template<typename U, typename C>
DistanceMatrix<T>* DistanceMatrixFactory<T>::getFilled(std::vector<U> *elements,
      const C &cmp) {
   DistanceMatrix<T>* dm = new DistanceMatrix<T> (elements->size());
     uint32_t numdist;
     uint32_t maxnumdist;

     numdist = 0;
     maxnumdist = (*dm).getNumElements() * (*dm).getNumElements();
     maxnumdist += (*dm).getNumElements();
     maxnumdist /= 2;
     //# pragma omp single
     for (uint32_t i = 0; i < (*dm).getNumElements(); i++) {
        #pragma omp parallel for shared(numdist)
        for (uint32_t j = i + 1; j < (*dm).getNumElements(); j++) {
           (*dm)(i, j) = cmp((*elements)[i], (*elements)[j]);

           # pragma omp atomic
           numdist++;

           if (numdist % 1000 == 0) {
              Log::dbg("distance matrix calculated %u/%u distances", numdist,
                    maxnumdist);
           }
        }
     }
     return dm;
}

#endif /*DISTMATRIX_CPP_*/

