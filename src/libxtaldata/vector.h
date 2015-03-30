#ifndef VECTOR_H_
#define VECTOR_H_

#include <string>
#include <vector>
#include <iterator>

#ifdef __SSE4_1__
#define _SSEVER __SSE4_1__
#endif

#ifdef _SSEVER
#include <smmintrin.h>
#endif

/*----------------------------------------------------------------------------*/

class Vector {
public:

   Vector();
   Vector(float x, float  y, float z);
   virtual ~Vector();

   float & operator[](int i);
   float operator[](int i) const;
   float* v();
   float & x();
   float & y();
   float & z();

   const float* v() const;
   float x() const;
   float y() const;
   float z() const;


   float dist(const Vector &c2) const;
   float length() const;
   float squared() const;
   void normalise();
   Vector normalised() const;
   float squareDist(const Vector &c2) const;
   Vector & operator+=(const Vector &c2);
   Vector & operator-=(const Vector &c2);
   Vector operator+(const Vector &c2) const;
   Vector operator-(const Vector &c2) const;

   Vector & operator+=(const float scale);
   Vector & operator-=(const float scale);
   Vector & operator*=(const float scale);
   Vector & operator/=(const float scale);
   Vector operator+(float scale) const;
   Vector operator-(float scale) const;
   Vector operator*(float scale) const;
   Vector operator/(float scale) const;

   virtual bool operator==(const Vector &v) const;
   virtual bool operator!=(const Vector &v) const;
   bool operator<(const Vector &v) const;


   Vector orthAxis() const;
   static Vector crossProd(const Vector &a, const Vector &b);
   static float dotProd(const Vector &a, const Vector &b);
   static float angle(const Vector &a, const Vector &b);
   static float distToPlane(const Vector &normal, const Vector &pointplane,
         const Vector &point);

   std::string toString() const;

   static Vector centre(const std::vector<Vector> &points);

protected:

   typedef union {
      float  f[4];
#ifdef _SSEVER
      __v4sf p;
#endif
   } VType;

   VType m_v;
};

#endif /* VECTOR_H_ */
