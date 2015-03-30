#ifndef DECOYREADERZD_H_
#define DECOYREADERZD_H_

#include "decoyreader.h"
#include <fstream>

/*----------------------------------------------------------------------------*/

class DecoyReaderZdPriv;

/*----------------------------------------------------------------------------*/

class DecoyReaderZd : public DecoyReader {
public:

   DecoyReaderZd(Mol *rec, Mol *lig, std::string zdFilename);
   virtual ~DecoyReaderZd();
   bool readNext(Mol *mol);
   Mol* getMovingMol();
   Mol* getFixedMol();
   const Matrix<float>* getMovRot();
   const Vector* getMovTrans();
   const Matrix<float>* getLigRot();
   const Vector* getLigTrans();
   bool good();
   uint32_t getNumDecoys();

private:
   void parseLine(Mol* mol, std::string line);
   DecoyReaderZdPriv *m;
};


#endif /* DECOYREADERZD_H_ */
