#include "kabschwrapper.h"
#include "kabsch.h"
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------*/

template<class T>
static void setMatrix(const T &a, const T &b, gsl_matrix *r, gsl_vector *t) {
   /* validate */
   if ((a.size() == 0) || (b.size() == 0)) {
      return;
   }
   /* setup matrices a and b */
   gsl_matrix *ma = gsl_matrix_alloc(a.size(), 3);
   for (unsigned int i = 0; i < a.size(); i++) {
      gsl_matrix_set(ma, i, 0, a[i]->x());
      gsl_matrix_set(ma, i, 1, a[i]->y());
      gsl_matrix_set(ma, i, 2, a[i]->z());
   }
   gsl_matrix *mb = gsl_matrix_alloc(b.size(), 3);
   for (unsigned int i = 0; i < b.size(); i++) {
      gsl_matrix_set(mb, i, 0, b[i]->x());
      gsl_matrix_set(mb, i, 1, b[i]->y());
      gsl_matrix_set(mb, i, 2, b[i]->z());
   }
   /* apply */
   kabsch(b.size(), mb, ma, r, t, NULL);
   gsl_matrix_free(ma);
   gsl_matrix_free(mb);
}

/*----------------------------------------------------------------------------*/

static void moveVectors(std::vector<Vector> &v, gsl_matrix *r, gsl_vector *t) {
   double x, y, z;
   for (unsigned int i = 0; i < v.size(); i++) {
      x = v[i].x() * gsl_matrix_get(r, 0, 0) + v[i].y() * gsl_matrix_get(r,
            0, 1) + v[i].z() * gsl_matrix_get(r, 0, 2) + gsl_vector_get(t, 0);
      y = v[i].x() * gsl_matrix_get(r, 1, 0) + v[i].y() * gsl_matrix_get(r,
            1, 1) + v[i].z() * gsl_matrix_get(r, 1, 2) + gsl_vector_get(t, 1);
      z = v[i].x() * gsl_matrix_get(r, 2, 0) + v[i].y() * gsl_matrix_get(r,
            2, 1) + v[i].z() * gsl_matrix_get(r, 2, 2) + gsl_vector_get(t, 2);
      v[i].x() = x;
      v[i].y() = y;
      v[i].z() = z;
   }
}

/*----------------------------------------------------------------------------*/
/* makes a copy of array b to calc rmsd */
float KabschWrapper::getCaRmsd(const Mol &a, const Mol &b) {
   std::vector<Vector> buffer;
   float rmsd;

   assert(a.size() == b.size());
   /* do kabsch */
   gsl_matrix *r = gsl_matrix_alloc(3, 3);
   gsl_vector *t = gsl_vector_alloc(3);
   setMatrix(a.atoms(), b.atoms(), r, t);
   //   /* make copy of copy bs vectors */
   buffer.resize(b.size());
   for (unsigned int i = 0; i < buffer.size(); i++) {
      buffer[i] = *b.atoms()[i];
   }
   /* move vectors to new position */
   moveVectors(buffer, r, t);
   /* free matrix and vector */
   gsl_matrix_free(r);
   gsl_vector_free(t);
   /* get rmsd */
   rmsd = 0;
   for (unsigned int i = 0; i < buffer.size(); i++) {
      if (a.atoms()[i]->isCa() == true) {
         rmsd += a.atoms()[i]->squareDist(buffer[i]);
      }
   }
   return sqrt(rmsd / buffer.size());
}

/*----------------------------------------------------------------------------*/

float KabschWrapper::getRmsd(const std::vector<Atom*> &a, const std::vector<
      Atom*> &b) {
	std::vector<Vector> buffer;
   float rmsd;

   assert(a.size() == b.size());
   /* do kabsch */
   gsl_matrix *r = gsl_matrix_alloc(3, 3);
   gsl_vector *t = gsl_vector_alloc(3);
   setMatrix(a, b, r, t);
   /* make copy of copy bs vectors */
   buffer.resize(b.size());
   for (unsigned int i = 0; i < buffer.size(); i++) {
      buffer[i] = *b[i];
   }
   /* move vectors to new position */
   moveVectors(buffer, r, t);
   /* free matrix and vector */
   gsl_matrix_free(r);
   gsl_vector_free(t);
   /* get rmsd */
   rmsd = 0;
   for (unsigned int i = 0; i < buffer.size(); i++) {
      if (a[i]->isCa() == true) {
         rmsd += a[i]->squareDist(buffer[i]);
      }
   }
   return sqrt(rmsd / buffer.size());
}

/*----------------------------------------------------------------------------*/

static void moveAtoms(std::vector<Atom *> &a, gsl_matrix *r, gsl_vector *t) {
   double x, y, z;
   for (unsigned int i = 0; i < a.size(); i++) {
      x = a[i]->x() * gsl_matrix_get(r, 0, 0) + a[i]->y()
            * gsl_matrix_get(r, 0, 1) + a[i]->z()
            * gsl_matrix_get(r, 0, 2) + gsl_vector_get(t, 0);
      y = a[i]->x() * gsl_matrix_get(r, 1, 0) + a[i]->y()
            * gsl_matrix_get(r, 1, 1) + a[i]->z()
            * gsl_matrix_get(r, 1, 2) + gsl_vector_get(t, 1);
      z = a[i]->x() * gsl_matrix_get(r, 2, 0) + a[i]->y()
            * gsl_matrix_get(r, 2, 1) + a[i]->z()
            * gsl_matrix_get(r, 2, 2) + gsl_vector_get(t, 2);
      a[i]->x() = x;
      a[i]->y() = y;
      a[i]->z() = z;
   }
}

/*----------------------------------------------------------------------------*/

void KabschWrapper::superimpose(std::vector<Atom*> &atoms, const std::vector<
      Atom*> &from, const std::vector<Atom*> &to, std::vector<float> *trans) {
   assert(from.size() == to.size());
   /* do kabsch */
   gsl_matrix *r = gsl_matrix_alloc(3, 3);
   gsl_vector *t = gsl_vector_alloc(3);
   setMatrix(to, from, r, t);
   moveAtoms(atoms, r, t);

   if (trans != NULL) {
      trans->clear();
      trans->push_back(gsl_matrix_get(r, 0, 0));
      trans->push_back(gsl_matrix_get(r, 0, 1));
      trans->push_back(gsl_matrix_get(r, 0, 2));
      trans->push_back(gsl_matrix_get(r, 1, 0));
      trans->push_back(gsl_matrix_get(r, 1, 1));
      trans->push_back(gsl_matrix_get(r, 1, 2));
      trans->push_back(gsl_matrix_get(r, 2, 0));
      trans->push_back(gsl_matrix_get(r, 2, 1));
      trans->push_back(gsl_matrix_get(r, 2, 2));
      trans->push_back(gsl_vector_get(t, 0));
      trans->push_back(gsl_vector_get(t, 1));
      trans->push_back(gsl_vector_get(t, 2));
   }

   /* free matrix and vector */
   gsl_matrix_free(r);
   gsl_vector_free(t);
}
