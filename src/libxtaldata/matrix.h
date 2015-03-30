#ifndef MATRIX_H_
#define MATRIX_H_

#include "vector.h"
#include <vector>
/*----------------------------------------------------------------------------*/

template <class T>
class Matrix : public std::vector<std::vector<T> > {
private:
   unsigned int m;
   unsigned int n;
public:
   Matrix(unsigned int m, unsigned int n, T init = 0);
   //void print() const;

   unsigned int getM() const;
   unsigned int getN() const;


//   bool write(std::ofstream *outfile) const;
//   bool write(std::string filename) const;
//
//   bool read(std::ifstream infile);
//   bool read(std::string filename);

   /** @brief multiplies a matrix with a vector */
   Vector const operator*(const Vector & c) const;

   Matrix<T> & operator*=(const float s);

   Matrix<T> & operator=(const Matrix<T> &m2);

   bool operator==(const Matrix<float> &m2) const;

   /* rot = r2 * r1; r1 is first r2 is second rotation */
   Matrix<T> operator*(const Matrix<T> &m2) const;

   /** @brief generates a 3x3 rotation matrix, which rotates around the axis v
     * @param v unit vector representing the axis to rotate around
     * @param a angle of rotation
     * @return rotation matrix */
   static Matrix<T> getRotationMatrix(const Vector &v, T a);

   static Matrix<T> getEulerRotation(T a, T b, T c);

   void transpose();
};

#endif /*MATRIX_H_*/
