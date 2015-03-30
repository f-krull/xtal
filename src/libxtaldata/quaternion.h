#ifndef QUATERNION_H_
#define QUATERNION_H_

#include "vector.h"
#include "matrix.cpp"

/*----------------------------------------------------------------------------*/

class Quaternion {
private:
protected:
public:

   Quaternion();
   Quaternion(float x, float y, float z, float w);
   Quaternion(const Matrix<float> &m);

   float & x();
   float & y();
   float & z();
   float & w();

   float x() const;
   float y() const;
   float z() const;
   float w() const;

   float *v();
   const float *v() const;

   void normalise();
   Quaternion getConjugate() const;
   Quaternion operator*(const Quaternion &rq) const;
   Quaternion & operator*=(const Quaternion &rq);
   Vector operator*(const Vector &vec) const;
   Matrix<float> getMatrix() const;
   void getMatrix(Matrix<float> &m) const;
   std::string toString() const;
   void getMatrix4v(float *matrix) const;
   static float innerProd(const Quaternion &q1, const Quaternion &q2);
   static Quaternion fromAxis(const Vector &v, float angle);
   static Quaternion fromEuler(float pitch, float yaw, float roll);
   float getRotAngle() const;

protected:
   float m_v[4];
};

#endif /*QUATERNION_H_*/
