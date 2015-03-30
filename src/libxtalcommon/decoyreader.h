#ifndef DECOYREADER_H_
#define DECOYREADER_H_

#include "../libxtaldata/mol.h"
#include "../libxtaldata/matrix.h"
#include <stdint.h>

/*----------------------------------------------------------------------------*/

class DecoyReader {
private:
protected:
public:
   virtual ~DecoyReader() {};
   virtual Mol* getMovingMol() = 0;
   virtual Mol* getFixedMol() = 0;
   virtual bool readNext(Mol *mov) = 0;
   virtual const Matrix<float>* getMovRot() = 0; /* move efficiently (smaller mol) */
   virtual const Vector* getMovTrans() = 0;
   virtual const Matrix<float>* getLigRot() = 0; /* move ligand */
   virtual const Vector* getLigTrans() = 0;
   virtual bool good() = 0;
   virtual uint32_t getNumDecoys() = 0;

   static DecoyReader* getDecoyReader(Mol *rec, Mol *lig, const char* filname);
};

#endif /* DECOYREADER_H_ */
