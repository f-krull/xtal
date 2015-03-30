#include "decoyreaderzd.h"
#include "../libxtalutil/log.h"
#include "../libxtaldata/matrix.cpp"
#include <assert.h>
#include <algorithm>

using namespace std;

/*----------------------------------------------------------------------------*/

class DecoyReaderZdPriv {
public:
   Mol* receptor;
   Mol* ligand;

   std::ifstream infile;
   std::string line;
   std::vector<Vector> pos_backup;
   Vector t, l, r, a, rand;
   float spacing, n;

   Matrix<float> m;
   Matrix<float> mrand;
   Matrix<float> rot;
   Vector trans;
   uint32_t numDecoys;

   DecoyReaderZdPriv();
};

/*----------------------------------------------------------------------------*/

DecoyReaderZdPriv::DecoyReaderZdPriv() :
   m(3, 3), mrand(3, 3), rot(3, 3) {
   numDecoys = 0;
}

/*----------------------------------------------------------------------------*/

DecoyReaderZd::DecoyReaderZd(Mol *rec, Mol *lig, string zdFilename) {
   m = new DecoyReaderZdPriv();

   m->receptor = rec;
   m->ligand = lig;
   /* make backup of coordinates */
   m->pos_backup.reserve(m->ligand->size());
   for (uint32_t i = 0; i < m->ligand->size(); i++) {
      m->pos_backup.push_back(*(*m->ligand).atoms()[i]);
   }
   /* open file for reading */
   m->infile.open(zdFilename.c_str());
   if (m->infile == NULL) {
      Log::err("ZdockReader cannot open file: %s", zdFilename.c_str());
      return;
   }
   string tmp;
   m->infile >> m->n;
   m->infile >> m->spacing;
   m->infile >> m->rand.x();
   m->infile >> m->rand.y();
   m->infile >> m->rand.z();
   m->infile >> tmp; /* rec */
   m->infile >> m->r.x();
   m->infile >> m->r.y();
   m->infile >> m->r.z();
   m->infile >> tmp; /* lig */
   m->infile >> m->l.x();
   m->infile >> m->l.y();
   m->infile >> m->l.z();

   m->mrand = Matrix<float>::getEulerRotation(m->rand.x(), m->rand.y(),
         m->rand.z());
   m->mrand.transpose();
}

/*----------------------------------------------------------------------------*/

DecoyReaderZd::~DecoyReaderZd() {
   m->infile.close();
   delete m;
}

/*----------------------------------------------------------------------------*/

bool DecoyReaderZd::good() {
   return m->infile.good();
}

/*----------------------------------------------------------------------------*/

static inline void moveMol(Mol &rotmol, vector<Vector> &pos_bup, const Matrix<
      float> &m, const Vector &t) {
   for (unsigned int i = 0; i < rotmol.size(); i++) {
      *rotmol.atoms()[i] = m * pos_bup[i];
      *rotmol.atoms()[i] += t;
   }
}

/*----------------------------------------------------------------------------*/

bool DecoyReaderZd::readNext(Mol *mol) {
   float ra;

   assert(mol->size() == m->pos_backup.size());
   while ((m->infile).fail() == false) {
      /* read 3 angles and 3 vectors */
      m->infile >> m->a.x();
      if ((m->infile).fail() == true) {
         return false;
      }
      m->infile >> m->a.y();
      m->infile >> m->a.z();
      m->infile >> m->t.x();
      m->infile >> m->t.y();
      m->infile >> m->t.z();
      m->infile >> ra; /* score */
      if ((m->infile).fail() == true) {
         return false;
      }
      /* rotation */
      m->m = Matrix<float>::getEulerRotation(m->a.x(), m->a.y(), m->a.z());
      m->m.transpose();
      if (m->t.x() >= m->n / 2) {
         m->t.x() -= m->n;
      }
      if (m->t.y() >= m->n / 2) {
         m->t.y() -= m->n;
      }
      if (m->t.z() >= m->n / 2) {
         m->t.z() -= m->n;
      }
      m->rot = m->m * m->mrand;
      m->trans = m->l;
      m->trans = m->rot * m->trans;
      m->trans *= -1;
      m->trans -= (m->t * m->spacing);
      m->trans += m->r;
      moveMol((*mol), m->pos_backup, m->rot, m->trans);
      return true;
   }
   return false;
}

/*----------------------------------------------------------------------------*/

Mol* DecoyReaderZd::getMovingMol() {
   return m->ligand;
}

/*----------------------------------------------------------------------------*/

Mol* DecoyReaderZd::getFixedMol() {
   return m->receptor;
}

/*----------------------------------------------------------------------------*/

const Matrix<float>* DecoyReaderZd::getMovRot() {
   return &m->rot;
}

/*----------------------------------------------------------------------------*/

const Vector* DecoyReaderZd::getMovTrans() {
   return &m->trans;
}

/*----------------------------------------------------------------------------*/

const Matrix<float>* DecoyReaderZd::getLigRot() {
   return &m->rot;
}

/*----------------------------------------------------------------------------*/

const Vector* DecoyReaderZd::getLigTrans() {
   return &m->trans;
}

/*----------------------------------------------------------------------------*/

uint32_t DecoyReaderZd::getNumDecoys() {
   if (m->numDecoys == 0) {
      istream::pos_type curr;
      curr = m->infile.tellg();
      m->infile.seekg(0, ios::beg);
      m->numDecoys = count(std::istreambuf_iterator<char>(m->infile),
            std::istreambuf_iterator<char>(), '\n');
      m->infile.seekg(curr);
      /* 4 lines header */
      if (m->numDecoys >= 4) {
         m->numDecoys -= 4;
      }
   }
   return m->numDecoys;
}

