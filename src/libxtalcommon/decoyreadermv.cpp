#include "decoyreadermv.h"
#include "../libxtalutil/log.h"
#include "../libxtaldata/matrix.cpp"
#include <assert.h>
#include <math.h>
#include <fstream>
#include <algorithm>

using namespace std;

/*----------------------------------------------------------------------------*/

class DecoyReaderMvPriv {
private:
protected:
public:
   DecoyReaderMvPriv();
   std::ifstream infile;
   int32_t numDecoys;

   Mol* receptor;
   Mol* ligand;
   Mol* bupmol;

   Matrix<float> rotMatrix;
   Vector trans;

};

/*----------------------------------------------------------------------------*/

DecoyReaderMvPriv::DecoyReaderMvPriv() :
   rotMatrix(3, 3) {
   numDecoys = 0;
}

/*----------------------------------------------------------------------------*/

DecoyReaderMv::DecoyReaderMv(Mol *rec, Mol *lig, const char* filename) {
   m = new DecoyReaderMvPriv();

   m->receptor = rec;
   m->ligand = lig;

   m->bupmol = new Mol(*m->ligand);
   m->infile.open(filename);
   if (!m->infile) {
      Log::err("DecoyReader cannot read decoy file: %s", filename);
      return;
   }
}

/*----------------------------------------------------------------------------*/

DecoyReaderMv::~DecoyReaderMv() {
   delete m->bupmol;
   m->infile.close();
   delete m;
}

/*----------------------------------------------------------------------------*/

static inline void readMatrix(ifstream *infile, Matrix<float> *m) {
   (*infile) >> (*m)[0][0];
   (*infile) >> (*m)[0][1];
   (*infile) >> (*m)[0][2];
   (*infile) >> (*m)[1][0];
   (*infile) >> (*m)[1][1];
   (*infile) >> (*m)[1][2];
   (*infile) >> (*m)[2][0];
   (*infile) >> (*m)[2][1];
   (*infile) >> (*m)[2][2];
}

/*----------------------------------------------------------------------------*/

static inline void readTrans(ifstream *infile, Vector *trans) {
   (*infile) >> (*trans).x();
   (*infile) >> (*trans).y();
   (*infile) >> (*trans).z();
}

/*----------------------------------------------------------------------------*/

static inline void moveMol(Mol &rotmol, Mol &bupmol, Matrix<float> *m,
      Vector *t) {
   for (unsigned int i = 0; i < rotmol.size(); i++) {
      *rotmol.atoms()[i] = (*m) * *bupmol.atoms()[i];
      *rotmol.atoms()[i] += (*t);
   }
}

/*----------------------------------------------------------------------------*/

bool DecoyReaderMv::readNext(Mol *mol) {
   assert(mol == m->ligand);
   readMatrix(&m->infile, &m->rotMatrix);
   readTrans(&m->infile, &m->trans);
   if (m->infile.fail()) {
      return false;
   }
   moveMol(*mol, *m->bupmol, &m->rotMatrix, &m->trans);
   return true;
}

/*----------------------------------------------------------------------------*/

Mol* DecoyReaderMv::getMovingMol() {
   return m->ligand;
}

/*----------------------------------------------------------------------------*/

Mol* DecoyReaderMv::getFixedMol() {
   return m->receptor;
}

/*----------------------------------------------------------------------------*/

const Matrix<float>* DecoyReaderMv::getMovRot() {
   return &m->rotMatrix;
}

/*----------------------------------------------------------------------------*/

const Vector* DecoyReaderMv::getMovTrans() {
   return &m->trans;
}

/*----------------------------------------------------------------------------*/

const Matrix<float>* DecoyReaderMv::getLigRot() {
   return &m->rotMatrix;
}

/*----------------------------------------------------------------------------*/

const Vector* DecoyReaderMv::getLigTrans() {
   return &m->trans;
}

/*----------------------------------------------------------------------------*/

bool DecoyReaderMv::good() {
   return m->infile.good();
}

/*----------------------------------------------------------------------------*/

uint32_t DecoyReaderMv::getNumDecoys() {
   if (m->numDecoys == 0) {
      istream::pos_type curr;
      curr = m->infile.tellg();
      m->infile.seekg(0, ios::beg);
      m->numDecoys = count(std::istreambuf_iterator<char>(m->infile),
            std::istreambuf_iterator<char>(), '\n');
      m->infile.seekg(curr);
   }
   return m->numDecoys;
}
