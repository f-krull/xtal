#ifndef DECOYREADERQ_H_
#define DECOYREADERQ_H_

#include "decoyreader.h"

class DecoyReaderQPriv;

/*----------------------------------------------------------------------------*/

class DecoyReaderQ : public DecoyReader {
public:
   DecoyReaderQ(Mol *rec, Mol *lig, const char* filename);
   virtual ~DecoyReaderQ();
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
   DecoyReaderQPriv *m;
};


#endif /* DECOYREADERQ_H_ */
