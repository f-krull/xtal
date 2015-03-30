#include "decoyreaderq.h"
#include "../libxtalutil/log.h"
#include "../libxtaldata/matrix.cpp"
#include "../libxtaldata/quaternion.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

/*----------------------------------------------------------------------------*/

typedef struct Header {
   char name[5];
   uint8_t version;
   uint8_t unused[250];
} Header;

/*----------------------------------------------------------------------------*/

//static inline void initHeader(Header &header) {
//   memcpy(header.name, "DECOY", 5);
//   header.version = 1;
//   memset(header.unused, 0, sizeof(header.unused));
//}

/*----------------------------------------------------------------------------*/

class DecoyReaderQPriv {
private:
protected:
public:
   DecoyReaderQPriv();
   FILE* infile;

   Mol* receptor;
   Mol* ligand;
   Mol* bupmol;

   Matrix<float> rotMatrix;
   Quaternion rotQ;
   Vector trans;
   uint32_t numDecoys;
};

/*----------------------------------------------------------------------------*/

inline bool writeDecoy(const Quaternion &q, const Vector &t, FILE *f) {
   static int8_t sep[4];
   bool err = false;

   memset(sep, 0, sizeof(sep));
   err = err || (fwrite(q.v(), sizeof(float_t), 4, f) != 4);
   err = err || (fwrite(t.v(), sizeof(float_t), 3, f) != 3);
   err = err || (fwrite(&sep, sizeof(uint8_t), 4, f) != 4);
   return !err;
}

/*----------------------------------------------------------------------------*/

inline bool readDecoy(Quaternion &q, Vector &t, FILE *f) {
   static int8_t sep[4];

   bool err = feof(f);
   err = err || (fread(q.v(), sizeof(float_t), 4, f) != 4);
   err = err || (fread(t.v(), sizeof(float_t), 3, f) != 3);
   err = err || (fread(&sep, sizeof(uint8_t), 4, f) != 4);
   return !err;
}

/*----------------------------------------------------------------------------*/

DecoyReaderQPriv::DecoyReaderQPriv() :
   rotMatrix(3, 3) {
   infile = NULL;
   numDecoys = 0;
}

/*----------------------------------------------------------------------------*/

DecoyReaderQ::DecoyReaderQ(Mol *rec, Mol *lig, const char* filename) {
   Header header;

   m = new DecoyReaderQPriv();

   m->receptor = rec;
   m->ligand = lig;

   m->bupmol = new Mol(*m->ligand);
   m->infile = fopen(filename, "r");
   if (m->infile == NULL) {
      Log::err("DecoyReader cannot read decoy file: %s", filename);
      return;
   }
   if (fread(&header, sizeof(header), 1, m->infile) != 1) {
      Log::err("DecoyReader cannot read decoy file: %s", filename);
      fclose(m->infile);
      m->infile = NULL;
      return;
   }
   if (strncmp(header.name, "DECOY", 5)) {
      Log::err("DecoyReader wrong file format %s\n", filename);
      fclose(m->infile);
      m->infile = NULL;
      return;
   }
}

/*----------------------------------------------------------------------------*/

DecoyReaderQ::~DecoyReaderQ() {
   delete m->bupmol;
   fclose(m->infile);
   delete m;
}

/*----------------------------------------------------------------------------*/

static inline void moveMol(Mol &rotmol, Mol &bupmol, const Matrix<float> &m,
      const Vector &t) {
   for (unsigned int i = 0; i < rotmol.size(); i++) {
      *rotmol.atoms()[i] = m * *bupmol.atoms()[i];
      *rotmol.atoms()[i] += t;
   }
}

/*----------------------------------------------------------------------------*/

bool DecoyReaderQ::readNext(Mol *mol) {
   assert(mol == m->ligand);
   if ((m->infile == NULL) || (readDecoy(m->rotQ, m->trans, m->infile) != true)) {
      return false;
   }
   m->rotQ.getMatrix(m->rotMatrix);
   moveMol(*mol, *m->bupmol, m->rotMatrix, m->trans);
   return true;
}

/*----------------------------------------------------------------------------*/

Mol* DecoyReaderQ::getMovingMol() {
   return m->ligand;
}

/*----------------------------------------------------------------------------*/

Mol* DecoyReaderQ::getFixedMol() {
   return m->receptor;
}

/*----------------------------------------------------------------------------*/

const Matrix<float>* DecoyReaderQ::getMovRot() {
   return &m->rotMatrix;
}

/*----------------------------------------------------------------------------*/

const Vector* DecoyReaderQ::getMovTrans() {
   return &m->trans;
}

/*----------------------------------------------------------------------------*/

const Matrix<float>* DecoyReaderQ::getLigRot() {
   return &m->rotMatrix;
}

/*----------------------------------------------------------------------------*/

const Vector* DecoyReaderQ::getLigTrans() {
   return &m->trans;
}

/*----------------------------------------------------------------------------*/

bool DecoyReaderQ::good() {
   return (m->infile != NULL) && (feof(m->infile) == 0);
}

/*----------------------------------------------------------------------------*/

uint32_t DecoyReaderQ::getNumDecoys() {
   if (m->numDecoys == 0) {
      long int curr, size;
      curr = ftell(m->infile);
      fseek(m->infile, 0L, SEEK_END);
      size = ftell(m->infile);
      fseek(m->infile, curr, SEEK_SET);
      m->numDecoys = (size - sizeof(Header)) / 32;
   }
   return m->numDecoys;
}
