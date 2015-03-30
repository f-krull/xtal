#include "./quaternion.h"
#include "../libxtalutil/common.h"
#include <cmath>
#include <sstream>
#include <memory.h>

/*----------------------------------------------------------------------------*/

Quaternion::Quaternion(float x, float y, float z, float w) {
   this->x() = x;
   this->y() = y;
   this->z() = z;
   this->w() = w;
}

/*----------------------------------------------------------------------------*/

Quaternion::Quaternion() {
   memset(m_v, 0, sizeof(m_v));
}

/*----------------------------------------------------------------------------*/

Quaternion::Quaternion(const Matrix<float> &m) {
   float t, s;

   t = 1 + m[0][0] + m[1][1] + m[2][2];

   if (t > 0.00000001) {
      s = sqrt(t) * 2;
      x() = (m[2][1] - m[1][2]) / s;
      y() = (m[0][2] - m[2][0]) / s;
      z() = (m[1][0] - m[0][1]) / s;
      w() = 0.25 * s;
   } else if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
      s = sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]) * 2;
      x() = 0.25 * s;
      y() = (m[0][1] + m[1][0]) / s;
      z() = (m[0][2] + m[2][0]) / s;
      w() = (m[2][1] - m[1][2]) / s;
   } else if (m[1][1] > m[2][2]) {
      s = sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]) * 2;
      x() = (m[1][0] + m[0][1]) / s;
      y() = 0.25 * s;
      z() = (m[2][1] + m[1][2]) / s;
      w() = (m[0][2] - m[2][0]) / s;
   } else {
      s = sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]) * 2;
      x() = (m[0][2] + m[2][0]) / s;
      y() = (m[2][1] + m[1][2]) / s;
      z() = 0.25 * s;
      w() = (m[1][0] - m[0][1]) / s;
   }

   //   X = ( mat[9] - mat[6] ) / s;
   //   Y = ( mat[2] - mat[8] ) / s;
   //   Z = ( mat[4] - mat[1] ) / s;

   //   0 1 2 3
   //   4 5 6 7
   //   8 9 0 1
   //   2 3 4 5

}

/*----------------------------------------------------------------------------*/

float & Quaternion::x() {
   return m_v[0];
}

/*----------------------------------------------------------------------------*/

float & Quaternion::y() {
   return m_v[1];
}

/*----------------------------------------------------------------------------*/


float & Quaternion::z() {
   return m_v[2];
}

/*----------------------------------------------------------------------------*/

float & Quaternion::w() {
   return m_v[3];
}

/*----------------------------------------------------------------------------*/

float Quaternion::x() const {
   return m_v[0];
}
/*----------------------------------------------------------------------------*/

float Quaternion::y() const {
   return m_v[1];
}

/*----------------------------------------------------------------------------*/

float Quaternion::z() const {
   return m_v[2];
}

/*----------------------------------------------------------------------------*/

float Quaternion::w() const {
   return m_v[3];
}

/*----------------------------------------------------------------------------*/

float * Quaternion::v() {
   return m_v;
}

/*----------------------------------------------------------------------------*/

const float * Quaternion::v() const {
   return m_v;
}

/*----------------------------------------------------------------------------*/

Quaternion Quaternion::getConjugate() const {
   return Quaternion(-x(), -y(), -z(), w());
}

/*----------------------------------------------------------------------------*/

Quaternion Quaternion::operator*(const Quaternion &rq) const {
   return Quaternion(w() * rq.x() + x() * rq.w() + y() * rq.z() - z() * rq.y(), w() * rq.y() + y()
         * rq.w() + z() * rq.x() - x() * rq.z(), w() * rq.z() + z() * rq.w() + x() * rq.y() - y()
         * rq.x(), w() * rq.w() - x() * rq.x() - y() * rq.y() - z() * rq.z());
}

/*----------------------------------------------------------------------------*/
Quaternion & Quaternion::operator*=(const Quaternion &rq) {
   float tmpx, tmpy, tmpz, tmpw;

   tmpx = w() * rq.x() + x() * rq.w() + y() * rq.z() - z() * rq.y();
   tmpy = w() * rq.y() + y() * rq.w() + z() * rq.x() - x() * rq.z();
   tmpz = w() * rq.z() + z() * rq.w() + x() * rq.y() - y() * rq.x();
   tmpw = w() * rq.w() - x() * rq.x() - y() * rq.y() - z() * rq.z();
   this->x() = tmpx;
   this->y() = tmpy;
   this->z() = tmpz;
   this->w() = tmpw;
   return (*this);
}

/*----------------------------------------------------------------------------*/

Vector Quaternion::operator*(const Vector &vec) const {
   Quaternion vecQuat, resQuat;
   Vector vn(vec);

   vn.normalise();
   vecQuat.x() = vn.x();
   vecQuat.y() = vn.y();
   vecQuat.z() = vn.z();
   vecQuat.w() = 0.0f;
   resQuat = vecQuat * getConjugate();
   resQuat = *this * resQuat;
   return (Vector(resQuat.x(), resQuat.y(), resQuat.z()));
}

/*----------------------------------------------------------------------------*/

Quaternion Quaternion::fromAxis(const Vector &v, float angle) {
   Quaternion q;
   float sinAngle;
   Vector vn(v);

   angle *= 0.5f;
   vn.normalise();
   sinAngle = sin(angle);
   q.x() = (vn.x() * sinAngle);
   q.y() = (vn.y() * sinAngle);
   q.z() = (vn.z() * sinAngle);
   q.w() = cos(angle);
   return q;
}

/*----------------------------------------------------------------------------*/

float Quaternion::getRotAngle() const {
   return acos(w()) * 2;
}

/*----------------------------------------------------------------------------*/

void Quaternion::normalise() {
   float mag;

   mag = w() * w() + x() * x() + y() * y() + z() * z();
   if (fabs(mag - 1.0f) > 0.000001f) {
      mag = sqrt(mag);
      w() /= mag;
      x() /= mag;
      y() /= mag;
      z() /= mag;
   }
}

/*----------------------------------------------------------------------------*/

Quaternion Quaternion::fromEuler(float pitch, float yaw, float roll) {
   float p, y, r, sinp, siny, sinr, cosp, cosy, cosr;
   Quaternion q;

   p = pitch * common::PI_180TH / 2.0;
   y = yaw * common::PI_180TH / 2.0;
   r = roll * common::PI_180TH / 2.0;
   sinp = sin(p);
   siny = sin(y);
   sinr = sin(r);
   cosp = cos(p);
   cosy = cos(y);
   cosr = cos(r);
   q.x() = sinr * cosp * cosy - cosr * sinp * siny;
   q.y() = cosr * sinp * cosy + sinr * cosp * siny;
   q.z() = cosr * cosp * siny - sinr * sinp * cosy;
   q.w() = cosr * cosp * cosy + sinr * sinp * siny;
   q.normalise();
   return q;
}

/*----------------------------------------------------------------------------*/

Matrix<float> Quaternion::getMatrix() const {
   float x2, y2, z2, xy, xz, yz, wx, wy, wz;
   Matrix<float> m(3, 3);

   x2 = x() * x();
   y2 = y() * y();
   z2 = z() * z();
   xy = x() * y();
   xz = x() * z();
   yz = y() * z();
   wx = w() * x();
   wy = w() * y();
   wz = w() * z();
   //   m[0][0] = 1.0f - 2.0f * (y2 + z2);
   //   m[1][0] = 2.0f * (xy - wz);
   //   m[2][0] = 2.0f * (xz + wy);
   //   m[0][1] = 2.0f * (xy + wz);
   //   m[1][1] = 1.0f - 2.0f * (x2 + z2);
   //   m[2][1] = 2.0f * (yz - wx);
   //   m[0][2] = 2.0f * (xz - wy);
   //   m[1][2] = 2.0f * (yz + wx);
   //   m[2][2] = 1.0f - 2.0f * (x2 + y2);
   m[0][0] = 1.0f - 2.0f * (y2 + z2);
   m[0][1] = 2.0f * (xy - wz);
   m[0][2] = 2.0f * (xz + wy);
   m[1][0] = 2.0f * (xy + wz);
   m[1][1] = 1.0f - 2.0f * (x2 + z2);
   m[1][2] = 2.0f * (yz - wx);
   m[2][0] = 2.0f * (xz - wy);
   m[2][1] = 2.0f * (yz + wx);
   m[2][2] = 1.0f - 2.0f * (x2 + y2);
   return m;
}

void Quaternion::getMatrix(Matrix<float> &m) const {
   float x2, y2, z2, xy, xz, yz, wx, wy, wz;

   assert(m.getM() == 3 && m.getN() == 3);
   x2 = x() * x();
   y2 = y() * y();
   z2 = z() * z();
   xy = x() * y();
   xz = x() * z();
   yz = y() * z();
   wx = w() * x();
   wy = w() * y();
   wz = w() * z();
   //   m[0][0] = 1.0f - 2.0f * (y2 + z2);
   //   m[1][0] = 2.0f * (xy - wz);
   //   m[2][0] = 2.0f * (xz + wy);
   //   m[0][1] = 2.0f * (xy + wz);
   //   m[1][1] = 1.0f - 2.0f * (x2 + z2);
   //   m[2][1] = 2.0f * (yz - wx);
   //   m[0][2] = 2.0f * (xz - wy);
   //   m[1][2] = 2.0f * (yz + wx);
   //   m[2][2] = 1.0f - 2.0f * (x2 + y2);
   m[0][0] = 1.0f - 2.0f * (y2 + z2);
   m[0][1] = 2.0f * (xy - wz);
   m[0][2] = 2.0f * (xz + wy);
   m[1][0] = 2.0f * (xy + wz);
   m[1][1] = 1.0f - 2.0f * (x2 + z2);
   m[1][2] = 2.0f * (yz - wx);
   m[2][0] = 2.0f * (xz - wy);
   m[2][1] = 2.0f * (yz + wx);
   m[2][2] = 1.0f - 2.0f * (x2 + y2);
}

/*----------------------------------------------------------------------------*/

void Quaternion::getMatrix4v(float *m) const {
   float x2, y2, z2, xy, xz, yz, wx, wy, wz;

   x2 = x() * x();
   y2 = y() * y();
   z2 = z() * z();
   xy = x() * y();
   xz = x() * z();
   yz = y() * z();
   wx = w() * x();
   wy = w() * y();
   wz = w() * z();
   m[0] = 1.0f - 2.0f * (y2 + z2);
   m[1] = 2.0f * (xy - wz);
   m[2] = 2.0f * (xz + wy);
   m[3] = 0.0f;
   m[4] = 2.0f * (xy + wz);
   m[5] = 1.0f - 2.0f * (x2 + z2);
   m[6] = 2.0f * (yz - wx);
   m[7] = 0.0f;
   m[8] = 2.0f * (xz - wy);
   m[9] = 2.0f * (yz + wx);
   m[10] = 1.0f - 2.0f * (x2 + y2);
   m[11] = 0.0f;
   m[12] = 0.0f;
   m[13] = 0.0f;
   m[14] = 0.0f;
   m[15] = 1.0f;
}

/*----------------------------------------------------------------------------*/

float Quaternion::innerProd(const Quaternion &q1, const Quaternion &q2) {
   return q1.x() * q2.x() + q1.y() * q2.y() + q1.z() * q2.z() + q1.w() * q2.w();
}

/*----------------------------------------------------------------------------*/

std::string Quaternion::toString() const {
   std::ostringstream ostr;

   ostr << x() << " " << y() << " " << z() << " " << w();
   return ostr.str();
}
