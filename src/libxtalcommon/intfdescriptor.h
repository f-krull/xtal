#ifndef INTFDESCRIPTOR_H_
#define INTFDESCRIPTOR_H_

#include <string>
#include <vector>
#include <stdint.h>
#include "../libxtaldata/protein.h"
#include "../libxtalcommon/interface.h"
typedef char cbool;


/*----------------------------------------------------------------------------*/

/* sequence containing interface info */
class UqChain {
public:
   /* sequence (upper case)*/
   std::string seq;
   /* interface (true) */
   std::vector<cbool> intf;
   /* interface size */
   uint32_t intfSize;

   std::string intfName;

   char name;

   UqChain(char name, const std::string &seq, const std::string &intfname);

   std::string getIntStr() const;

};

/*----------------------------------------------------------------------------*/

class UqEntry {
public:
   std::string intname;
   std::string name;
   float resolution;
   uint32_t year;
   std::string chains1;
   std::string chains2;
   std::vector<UqChain> seq1;
   std::vector<UqChain> seq2;
   std::string title;
   std::string ligands;
   uint32_t numNonHohLig;

   UqEntry();

   /* favor small resolutions, greater years and greater names */
   bool operator<(const UqEntry &e) const;

   bool read(const std::string &fnpBC, const std::string &cnB1,
         const std::string &cnB2, const IntfDef &id);


   bool read(const Protein &pBC, const Chains &cB1, const Chains &cB2, const IntfDef &id);

   void print();
};

/*----------------------------------------------------------------------------*/

class DistUqentryCmp {
public:

   DistUqentryCmp(bool showAlignment) : m_showAlignment(showAlignment) {
      m_ignoreSizeDiff = false;
   }

   float operator()(const UqEntry &e1, const UqEntry &e2) const;

   void setIgnoreSizeDiff(bool b) {m_ignoreSizeDiff = b;}

private:
   bool m_showAlignment;
   bool m_ignoreSizeDiff;
};


#endif /* INTFDESCRIPTOR_H_ */
