#include "vector.h"
#include <math.h>
#include <memory.h>
#include <assert.h>
#include <sstream>

/*----------------------------------------------------------------------------*/

Vector::Vector() {
   memset(m_v.f, 0, sizeof(m_v));
}

/*----------------------------------------------------------------------------*/

Vector::Vector(float x, float  y, float z) {
   m_v.f[0] = x;
   m_v.f[1] = y;
   m_v.f[2] = z;
}

/*----------------------------------------------------------------------------*/

Vector::~Vector() {
}

/*----------------------------------------------------------------------------*/

float* Vector::v() {
   return m_v.f;
}

/*----------------------------------------------------------------------------*/

float & Vector::x() {
   return m_v.f[0];
}

/*----------------------------------------------------------------------------*/

float & Vector::y() {
   return m_v.f[1];
}

/*----------------------------------------------------------------------------*/

float & Vector::z() {
   return m_v.f[2];
}

/*----------------------------------------------------------------------------*/

float & Vector::operator[](int i) {
   assert(i > -1 && i < 3);
   return m_v.f[i];
}

/*----------------------------------------------------------------------------*/

float Vector::operator[](int i) const {
   assert(i > -1 && i < 3);
   return m_v.f[i];
}

/*----------------------------------------------------------------------------*/

const float* Vector::v() const {
   return m_v.f;
}

/*----------------------------------------------------------------------------*/

float Vector::x() const{
   return m_v.f[0];
}

/*----------------------------------------------------------------------------*/

float Vector::y() const{
   return m_v.f[1];
}

/*----------------------------------------------------------------------------*/

float Vector::z() const{
   return m_v.f[2];
}

/*----------------------------------------------------------------------------*/

std::string Vector::toString() const {
   std::ostringstream ostr;

   ostr << m_v.f[0] << " " << m_v.f[1] << " " << m_v.f[2];
   return ostr.str();
}

/*----------------------------------------------------------------------------*/

float Vector::dist(const Vector &c2) const {
   return sqrt(this->squareDist(c2));
}

/*----------------------------------------------------------------------------*/

float Vector::length() const {
#ifndef _SSEVER
   return sqrt(pow(m_v.f[0], 2) + pow(m_v.f[1], 2) + pow(m_v.f[2], 2));
#else
//   VType t;
//   t.p = m_v.p * m_v.p;
//   return sqrt(t.f[0] + t.f[1] + t.f[2]);
   return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(m_v.p, m_v.p, 0x71)));
#endif
}

/*----------------------------------------------------------------------------*/

float Vector::squared() const {
#ifndef _SSEVER
   return pow(m_v.f[0], 2) + pow(m_v.f[1], 2) + pow(m_v.f[2], 2);
#else
//   VType t;
//   t.p = m_v.p * m_v.p;
//   return t.f[0] + t.f[1] + t.f[2];
   return _mm_cvtss_f32(_mm_dp_ps(m_v.p, m_v.p, 0x71));
#endif
}

/*----------------------------------------------------------------------------*/

void Vector::normalise() {
   (*this) /= length();
}

/*----------------------------------------------------------------------------*/

Vector Vector::normalised() const {
   return Vector(*this) /= length();
}

/*----------------------------------------------------------------------------*/

float Vector::squareDist(const Vector &c2) const {
#ifndef _SSEVER
   return pow(m_v.f[0] - c2.m_v.f[0], 2) + pow(m_v.f[1] - c2.m_v.f[1], 2) + pow(m_v.f[2]
         - c2.m_v.f[2], 2);
#else
   VType t;
   t.p = m_v.p - c2.m_v.p;
   return _mm_cvtss_f32(_mm_dp_ps(t.p, t.p, 0x71));
#endif
}

/*----------------------------------------------------------------------------*/

Vector & Vector::operator+=(const Vector &c2) {
#ifndef _SSEVER
   this->m_v.f[0] += c2.m_v.f[0];
   this->m_v.f[1] += c2.m_v.f[1];
   this->m_v.f[2] += c2.m_v.f[2];
#else
   m_v.p += c2.m_v.p;
#endif
   return (*this);
}

/*----------------------------------------------------------------------------*/

Vector & Vector::operator-=(const Vector &c2) {
#ifndef _SSEVER
   this->m_v.f[0] -= c2.m_v.f[0];
   this->m_v.f[1] -= c2.m_v.f[1];
   this->m_v.f[2] -= c2.m_v.f[2];
#else
   m_v.p -= c2.m_v.p;
#endif
   return (*this);
}

/*----------------------------------------------------------------------------*/

Vector Vector::operator+(const Vector &c2) const {
   return Vector(*this) += c2;
}

/*----------------------------------------------------------------------------*/

Vector Vector::operator-(const Vector &c2) const {
   return Vector(*this) -= c2;
}

/*----------------------------------------------------------------------------*/

Vector & Vector::operator+=(const float scale) {
   this->m_v.f[0] += scale;
   this->m_v.f[1] += scale;
   this->m_v.f[2] += scale;
   return (*this);
}

/*----------------------------------------------------------------------------*/

Vector & Vector::operator-=(const float scale) {
   this->m_v.f[0] -= scale;
   this->m_v.f[1] -= scale;
   this->m_v.f[2] -= scale;
   return (*this);
}

/*----------------------------------------------------------------------------*/

Vector & Vector::operator*=(const float scale) {
   this->m_v.f[0] *= scale;
   this->m_v.f[1] *= scale;
   this->m_v.f[2] *= scale;
   return (*this);
}

/*----------------------------------------------------------------------------*/

Vector & Vector::operator/=(const float scale) {
   this->m_v.f[0] /= scale;
   this->m_v.f[1] /= scale;
   this->m_v.f[2] /= scale;
   return (*this);
}

/*----------------------------------------------------------------------------*/

Vector Vector::operator+(float scale) const {
   Vector ret(*this);
   ret.m_v.f[0] += scale;
   ret.m_v.f[1] += scale;
   ret.m_v.f[2] += scale;
   return ret;
}

/*----------------------------------------------------------------------------*/

Vector Vector::operator-(float scale) const {
   Vector ret(*this);
   ret.m_v.f[0] -= scale;
   ret.m_v.f[1] -= scale;
   ret.m_v.f[2] -= scale;
   return ret;
}

/*----------------------------------------------------------------------------*/

Vector Vector::operator*(float scale) const {
   Vector ret(*this);
   ret.m_v.f[0] *= scale;
   ret.m_v.f[1] *= scale;
   ret.m_v.f[2] *= scale;
   return ret;
}

/*----------------------------------------------------------------------------*/

Vector Vector::operator/(float scale) const {
   Vector ret(*this);
   ret.m_v.f[0] /= scale;
   ret.m_v.f[1] /= scale;
   ret.m_v.f[2] /= scale;
   return ret;
}

/*----------------------------------------------------------------------------*/

Vector Vector::crossProd(const Vector &a, const Vector &b) {
   Vector cp;

   cp.m_v.f[0] = a.m_v.f[1] * b.m_v.f[2] - a.m_v.f[2] * b.m_v.f[1];
   cp.m_v.f[1] = a.m_v.f[2] * b.m_v.f[0] - a.m_v.f[0] * b.m_v.f[2];
   cp.m_v.f[2] = a.m_v.f[0] * b.m_v.f[1] - a.m_v.f[1] * b.m_v.f[0];
   return cp;
}

/*----------------------------------------------------------------------------*/

bool Vector::operator==(const Vector &v) const {
   if (fabs(this->m_v.f[0] - v.m_v.f[0]) > 0.0000001) {
      return false;
   }
   if (fabs(this->m_v.f[1] - v.m_v.f[1]) > 0.0000001) {
      return false;
   }
   if (fabs(this->m_v.f[2] - v.m_v.f[2]) > 0.0000001) {
      return false;
   }
   return true;
}

/*----------------------------------------------------------------------------*/

bool Vector::operator!=(const Vector &v) const {
   return !(*this == v);
}

/*----------------------------------------------------------------------------*/

bool Vector::operator<(const Vector &v) const {
   return this->length() < v.length();
}

/*----------------------------------------------------------------------------*/

float Vector::dotProd(const Vector &a, const Vector &b) {
   return (a.m_v.f[0] * b.m_v.f[0]) + (a.m_v.f[1] * b.m_v.f[1]) + (a.m_v.f[2] * b.m_v.f[2]);
}

/*----------------------------------------------------------------------------*/

float Vector::angle(const Vector &a, const Vector &b) {
   return acos(dotProd(a.normalised(), b.normalised()));
}

/*----------------------------------------------------------------------------*/

Vector Vector::orthAxis() const {
   Vector o;

   if ((m_v.f[2]) != 0) {
      o.m_v.f[0] = 1;
      o.m_v.f[1] = 1;
      o.m_v.f[2] = (m_v.f[0] * o.m_v.f[0] + m_v.f[1] * o.m_v.f[1]) * (-1) / m_v.f[2];
   } else if ((m_v.f[1]) != 0) {
      o.m_v.f[0] = 1;
      o.m_v.f[1] = (m_v.f[0] * o.m_v.f[0]) * (-1) / m_v.f[1];
      o.m_v.f[2] = 0;
   } else {
      /* only r.x != 0 */
      o.m_v.f[0] = 0;
      o.m_v.f[1] = 1;
      o.m_v.f[2] = 0;
   }
   return o;
}

/*----------------------------------------------------------------------------*/

float Vector::distToPlane(const Vector &normal, const Vector &pointplane,
      const Vector &point) {
   return Vector::dotProd(normal, point - pointplane);
}

/*----------------------------------------------------------------------------*/

Vector Vector::centre(const std::vector<Vector> &points) {
   Vector centre;
   for (unsigned int i = 0; i < points.size(); i++) {
      centre += points[i];
   }
   return (centre /= points.size());
}

/*----------------------------------------------------------------------------*/
//
//Vector::Vector(const Vector &v2) {
//   memcpy(m_v, v2.m_v, sizeof(m_v));
//}
//
/*----------------------------------------------------------------------------*/
//
//Vector & Vector::operator=(const Vector &v2) {
//   memcpy(m_v, v2.m_v, sizeof(m_v));
//   return (*this);
//}

