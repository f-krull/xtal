#ifndef DECOYREADERMV_H_
#define DECOYREADERMV_H_

#include "decoyreader.h"

class DecoyReaderMvPriv;

/*----------------------------------------------------------------------------*/

class DecoyReaderMv :public DecoyReader {
public:
   DecoyReaderMv(Mol *rec, Mol *lig, const char* filename);
   virtual ~DecoyReaderMv();
   bool readNext(Mol *lig);
   Mol* getMovingMol();
   Mol* getFixedMol();
   const Matrix<float>* getMovRot();
   const Vector* getMovTrans();
   const Matrix<float>* getLigRot();
   const Vector* getLigTrans();
   bool good();
   uint32_t getNumDecoys();
protected:
private:
   DecoyReaderMvPriv *m;
};

#endif /* DECOYREADERMV_H_ */
