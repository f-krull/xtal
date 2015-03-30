#ifndef AAMAP_H_
#define AAMAP_H_

#include <stdint.h>

/*----------------------------------------------------------------------------*/

class AaMapPriv;

/*----------------------------------------------------------------------------*/

class AaMap {
public:

   static unsigned char getId(char resName1);
   static unsigned char getId(const char* resName3);
   static bool isAa(const char* aa3str);
   static char getResn1(const char* resName3);
   static char idxToResn1(uint32_t idx);
   static const char* getResn3(unsigned char resId);
   ~AaMap();

   unsigned char _getId(char resName1) const;
   unsigned char _getId(const char* resName3) const;
   bool _isAa(const char* aa3str) const;
   char _getResn1(const char* resName3) const;
   char _idxToResn1(uint32_t idx) const;
   const char* _getResn3(unsigned char resId) const;


   static char UNK;
   static char ALA;
   static char CYS;
   static char ASP;
   static char GLU;
   static char PHE;
   static char GLY;
   static char HIS;
   static char ILE;
   static char LYS;
   static char LEU;
   static char MET;
   static char ASN;
   static char PRO;
   static char GLN;
   static char ARG;
   static char SER;
   static char THR;
   static char VAL;
   static char TRP;
   static char TYR;
   static char MSE;
   static char PYL;
   static char SEC;
   static char HSD;
   static char HSE;
   static char HSP;
   static char size;

private:
   void init();
   AaMap();
   AaMap(const AaMap &);
   AaMap & operator=(const AaMap &);
   static AaMap& getInstance();

   AaMapPriv *m;
};





#endif /*AAMAP_H_*/
