
#ifndef DISTMATRIX_H_
#define DISTMATRIX_H_

#include <vector>
#include <stdint.h>

/*----------------------------------------------------------------------------*/

template<typename T>
class DistanceMatrix {
public:
   DistanceMatrix(uint32_t numElements, T initValue = 0);
   ~DistanceMatrix();
   T get(uint32_t i, uint32_t j) const;
   void set(uint32_t i, uint32_t j, T value);
   void setAll(T value);
   uint32_t getNumElements();
   void print();
   T& operator()(uint32_t i, uint32_t j);

   bool save(const char* filename);
   bool load(const char* filename);
private:
   T ** m_matrix;
   T m_null;
protected:
   uint32_t m_numElements;
};

/*----------------------------------------------------------------------------*/

template <class T>
class DistanceMatrixFactory {
private:
protected:
public:
   template<typename U>
   static DistanceMatrix<T>* getFilled(std::vector<U> *elements,
         T(*dist)(const U&, const U&));

   template<typename U, typename C>
   static DistanceMatrix<T>* getFilled(std::vector<U> *elements,
         const C &cmp);
};

#endif /*DISTMATRIX_H_*/
