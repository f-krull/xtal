#ifndef MATRIX_CPP_
#define MATRIX_CPP_

#include "matrix.h"
#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------*/

template <class T>
Matrix<T>::Matrix(unsigned int w, unsigned int h, T init) {
   this->m = w;
   this->n = h;
   this->resize(w);
   for (unsigned int i = 0; i < this->size(); i++) {
      this->at(i).resize(h);
      for (unsigned int j = 0; j < this->at(i).size(); j++) {
         this->at(i)[j] = init;
      }
   }
}

/*----------------------------------------------------------------------------*/

//template <class T>
//void Matrix<T>::print() const {
//   for (unsigned int i = 0; i < this->size(); i++) {
//      for (unsigned int j = 0; j < this->at(i).size(); j++) {
//         std::cout << this->at(i)[j] << " ";
//      }
//      std::cout << std::endl;
//   }
//}

/*----------------------------------------------------------------------------*/

template <class T>
Matrix<T> Matrix<T>::getRotationMatrix(const Vector &v, T a) {
   Matrix m(3, 3);
   T cosa;
   T sina;

   sina = sin(a);
   cosa = cos(a);

   (m)[0][0] = cosa + (1 - cosa) * pow(v.x(), 2);
   (m)[0][1] = (1 - cosa) * v.x() * v.y() - sina * v.z();
   (m)[0][2] = (1 - cosa) * v.x() * v.z() + sina * v.y();

   (m)[1][0] = (1 - cosa) * v.y() * v.x() + sina * v.z();
   (m)[1][1] = cosa + (1 - cosa) * pow(v.y(), 2);
   (m)[1][2] = (1 - cosa) * v.y() * v.z() - sina * v.x();

   (m)[2][0] = (1 - cosa) * v.z() * v.x() - sina * v.y();
   (m)[2][1] = (1 - cosa) * v.z() * v.y() + sina * v.x();
   (m)[2][2] = cosa + (1 - cosa) * pow(v.z(), 2);
   return m;
}

/*----------------------------------------------------------------------------*/

template <class T>
Matrix<T> Matrix<T>::getEulerRotation(T phi, T the, T psi) {
   Matrix m(3, 3);
   T costhe, sinthe;
   T cosphi, sinphi;
   T cospsi, sinpsi;

   sinthe = sin(the);
   costhe = cos(the);
   sinpsi = sin(psi);
   cospsi = cos(psi);
   sinphi = sin(phi);
   cosphi = cos(phi);

   m[0][0] =  cospsi*cosphi - costhe*sinphi*sinpsi;
   m[0][1] =  cospsi*sinphi + costhe*cosphi*sinpsi;
   m[0][2] =  sinpsi*sinthe;
   m[1][0] = -sinpsi*cosphi - costhe*sinphi*cospsi;
   m[1][1] = -sinpsi*sinphi + costhe*cosphi*cospsi;
   m[1][2] =  cospsi*sinthe;
   m[2][0] =  sinthe*sinphi;
   m[2][1] = -sinthe*cosphi;
   m[2][2] =  costhe;

   return m;
}

/*----------------------------------------------------------------------------*/

template <class T>
Vector const Matrix<T>::operator*(const Vector & c) const {
   Vector c_tmp;

   /* TODO: not generic case yet */
   assert((*this).size() == 3);
   assert((*this)[0].size() == 3);
   c_tmp.x() = c.x() * (*this)[0][0]
             + c.y() * (*this)[0][1]
             + c.z() * (*this)[0][2];
   c_tmp.y() = c.x() * (*this)[1][0]
             + c.y() * (*this)[1][1]
             + c.z() * (*this)[1][2];
   c_tmp.z() = c.x() * (*this)[2][0]
             + c.y() * (*this)[2][1]
             + c.z() * (*this)[2][2];
   return c_tmp;
}

/*----------------------------------------------------------------------------*/

template <class T>
Matrix<T> & Matrix<T>::operator*=(const float s) {
   for (unsigned int i = 0; i < (*this).size(); i++) {
      for (unsigned int j = 0; j < (*this)[0].size(); j++) {
         (*this)[i][j] *= s;
      }
   }
   return (*this);
}

/*----------------------------------------------------------------------------*/

template <class T>
Matrix<T> & Matrix<T>::operator=(const Matrix<T> &m2) {
   assert((*this).size() == m2.size());
   assert((*this).at(0).size() == m2.at(0).size());
   for (unsigned int i = 0; i < (*this).size(); i++) {
      for (unsigned int j = 0; j < (*this)[0].size(); j++) {
         (*this)[i][j] = m2[i][j];
      }
   }
   return (*this);
}

/*----------------------------------------------------------------------------*/

template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &m2) const {
   Matrix<T> m(m2.size(), m2.at(0).size());

   /* only for 3x3 matrices so far */
   assert((*this).size() == 3);
   assert((*this).size() == m2.size());
   assert((*this).at(0).size() == m2.at(0).size());
   for (unsigned int i = 0; i < (*this).size(); i++) {
      for (unsigned int j = 0; j < (*this)[0].size(); j++) {
         m[i][j] = 0;
         for (unsigned int k = 0; k < 3; k++) {
            m[i][j] += (*this)[i][k] * m2[k][j];
         }
      }
   }
   return m;
}


/*----------------------------------------------------------------------------*/

template <class T>
void Matrix<T>::transpose() {
   T tmp;

   assert((*this).size() == (*this)[0].size());
   for (unsigned int i = 0; i < (*this).size(); i++) {
      for (unsigned int j = i + 1; j < (*this).size(); j++) {
         tmp = (*this)[i][j];
         (*this)[i][j] = (*this)[j][i];
         (*this)[j][i] = tmp;
      }
   }
}

/*----------------------------------------------------------------------------*/

template <class T>
unsigned int  Matrix<T>::getM() const {
   return m;
}

/*----------------------------------------------------------------------------*/

template <class T>
unsigned int  Matrix<T>::getN() const {
   return n;
}

/*----------------------------------------------------------------------------*/

template <class T>
bool Matrix<T>::operator==(const Matrix<float> &m2) const {
   if (this->getN() != m2.getN()) {
      return false;
   }
   if (this->getM() != m2.getM()) {
      return false;
   }
   for (unsigned int i = 0; i < getN(); i++) {
      for (unsigned int j = 0; j < getM(); j++) {
         if (fabs((*this)[i][j] - m2[i][j]) > 0.001) {
            return false;
         }
      }
   }
   return true;
}

/*----------------------------------------------------------------------------*/

//template <class T>
//bool Matrix<T>::write(std::ofstream *outfile) const {
//   (*outfile) << "matrix " << m << " " << n << "\n";
//   for (unsigned int i = 0; i < this->size(); i++) {
//      for (unsigned int j = 0; j < this->at(i).size(); j++) {
//         (*outfile) << this->at(i)[j] << " ";
//      }
//      (*outfile) << "\n";
//   }
//   return !((*outfile).fail());
//}
//
///*----------------------------------------------------------------------------*/
//
//template <class T>
//bool Matrix<T>::write(std::string outfilename) const {
//   std::ofstream outfile;
//
//   outfile.open(outfilename.c_str());
//   if (outfile == false) {
//      return false;
//   }
//   return write(&outfile);
//}

#endif /*MATRIX_CPP_*/
