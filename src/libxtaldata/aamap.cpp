#include "aamap.h"
#include "../libxtalutil/log.h"
#include <vector>
#include <map>
#include <string>
#include <assert.h>
#include <memory>

/*----------------------------------------------------------------------------*/

class AaMapPriv {
public:
   std::vector<char> aa1table;
   std::vector<std::string> aa3table;
   std::map<char, unsigned char> aa1map;
   std::map<std::string, unsigned char> aa3map;
   std::map<std::string, char> aa3toaa3map;
};

/*----------------------------------------------------------------------------*/

char AaMap::UNK = 0;
char AaMap::ALA = 1;
char AaMap::CYS = 2;
char AaMap::ASP = 3;
char AaMap::GLU = 4;
char AaMap::PHE = 5;
char AaMap::GLY = 6;
char AaMap::HIS = 7;
char AaMap::ILE = 8;
char AaMap::LYS = 9;
char AaMap::LEU = 10;
char AaMap::MET = 11;
char AaMap::ASN = 12;
char AaMap::PRO = 13;
char AaMap::GLN = 14;
char AaMap::ARG = 15;
char AaMap::SER = 16;
char AaMap::THR = 17;
char AaMap::VAL = 18;
char AaMap::TRP = 19;
char AaMap::TYR = 20;
char AaMap::size = 21;
/* alternatives: */
char AaMap::MSE = 11;
char AaMap::PYL = 10;
char AaMap::SEC = 0;
char AaMap::HSD = 7;
char AaMap::HSE = 7;
char AaMap::HSP = 7;

/*----------------------------------------------------------------------------*/

void AaMap::init() {

   m->aa1table.resize(21);
   m->aa3table.resize(21);

   m->aa3table[UNK] = "UNK";
   m->aa3table[ALA] = "ALA";
   m->aa3table[CYS] = "CYS";
   m->aa3table[ASP] = "ASP";
   m->aa3table[GLU] = "GLU";
   m->aa3table[PHE] = "PHE";
   m->aa3table[GLY] = "GLY";
   m->aa3table[HIS] = "HIS";
   m->aa3table[ILE] = "ILE";
   m->aa3table[LYS] = "LYS";
   m->aa3table[LEU] = "LEU";
   m->aa3table[MET] = "MET";
   m->aa3table[ASN] = "ASN";
   m->aa3table[PRO] = "PRO";
   m->aa3table[GLN] = "GLN";
   m->aa3table[ARG] = "ARG";
   m->aa3table[SER] = "SER";
   m->aa3table[THR] = "THR";
   m->aa3table[VAL] = "VAL";
   m->aa3table[TRP] = "TRP";
   m->aa3table[TYR] = "TYR";

   m->aa1table[UNK] = '?';
   m->aa1table[ALA] = 'A';
   m->aa1table[CYS] = 'C';
   m->aa1table[ASP] = 'D';
   m->aa1table[GLU] = 'E';
   m->aa1table[PHE] = 'F';
   m->aa1table[GLY] = 'G';
   m->aa1table[HIS] = 'H';
   m->aa1table[ILE] = 'I';
   m->aa1table[LYS] = 'K';
   m->aa1table[LEU] = 'L';
   m->aa1table[MET] = 'M';
   m->aa1table[ASN] = 'N';
   m->aa1table[PRO] = 'P';
   m->aa1table[GLN] = 'Q';
   m->aa1table[ARG] = 'R';
   m->aa1table[SER] = 'S';
   m->aa1table[THR] = 'T';
   m->aa1table[VAL] = 'V';
   m->aa1table[TRP] = 'W';
   m->aa1table[TYR] = 'Y';

   m->aa1map['?'] = UNK;
   m->aa1map['A'] = ALA;
   m->aa1map['C'] = CYS;
   m->aa1map['D'] = ASP;
   m->aa1map['E'] = GLU;
   m->aa1map['F'] = PHE;
   m->aa1map['G'] = GLY;
   m->aa1map['H'] = HIS;
   m->aa1map['I'] = ILE;
   m->aa1map['K'] = LYS;
   m->aa1map['L'] = LEU;
   m->aa1map['M'] = MET;
   m->aa1map['N'] = ASN;
   m->aa1map['P'] = PRO;
   m->aa1map['Q'] = GLN;
   m->aa1map['R'] = ARG;
   m->aa1map['S'] = SER;
   m->aa1map['T'] = THR;
   m->aa1map['V'] = VAL;
   m->aa1map['W'] = TRP;
   m->aa1map['Y'] = TYR;
   m->aa1map['U'] = SEC;

   m->aa3map["UNK"] = UNK;
   m->aa3map["ALA"] = ALA;
   m->aa3map["CYS"] = CYS;
   m->aa3map["ASP"] = ASP;
   m->aa3map["GLU"] = GLU;
   m->aa3map["PHE"] = PHE;
   m->aa3map["GLY"] = GLY;
   m->aa3map["HIS"] = HIS;
   m->aa3map["ILE"] = ILE;
   m->aa3map["LYS"] = LYS;
   m->aa3map["LEU"] = LEU;
   m->aa3map["MET"] = MET;
   m->aa3map["ASN"] = ASN;
   m->aa3map["PRO"] = PRO;
   m->aa3map["GLN"] = GLN;
   m->aa3map["ARG"] = ARG;
   m->aa3map["SER"] = SER;
   m->aa3map["THR"] = THR;
   m->aa3map["VAL"] = VAL;
   m->aa3map["TRP"] = TRP;
   m->aa3map["TYR"] = TYR;
   m->aa3map["MSE"] = MSE;
   m->aa3map["PYL"] = PYL;
   m->aa3map["SEC"] = SEC;
   m->aa3map["HSD"] = HSD;
   m->aa3map["HSE"] = HSE;
   m->aa3map["HSP"] = HSP;

   m->aa3toaa3map["UNK"] = '?';
   m->aa3toaa3map["ALA"] = 'A';
   m->aa3toaa3map["CYS"] = 'C';
   m->aa3toaa3map["ASP"] = 'D';
   m->aa3toaa3map["GLU"] = 'E';
   m->aa3toaa3map["PHE"] = 'F';
   m->aa3toaa3map["GLY"] = 'G';
   m->aa3toaa3map["HIS"] = 'H';
   m->aa3toaa3map["ILE"] = 'I';
   m->aa3toaa3map["LYS"] = 'K';
   m->aa3toaa3map["LEU"] = 'L';
   m->aa3toaa3map["MET"] = 'M';
   m->aa3toaa3map["ASN"] = 'N';
   m->aa3toaa3map["PRO"] = 'P';
   m->aa3toaa3map["GLN"] = 'Q';
   m->aa3toaa3map["ARG"] = 'R';
   m->aa3toaa3map["SER"] = 'S';
   m->aa3toaa3map["THR"] = 'T';
   m->aa3toaa3map["VAL"] = 'V';
   m->aa3toaa3map["TRP"] = 'W';
   m->aa3toaa3map["TYR"] = 'Y';
   m->aa3toaa3map["MSE"] = 'M';
   m->aa3toaa3map["PYL"] = 'L';
   m->aa3toaa3map["HSD"] = 'H';
   m->aa3toaa3map["HSE"] = 'H';
   m->aa3toaa3map["HSP"] = 'H';
}

/*----------------------------------------------------------------------------*/

unsigned char AaMap::_getId(char aa1) const {
   unsigned char res = AaMap::UNK;
   {
      std::map<char, unsigned char>::const_iterator it = m->aa1map.find(aa1);
      if (it != m->aa1map.end()) {
         res = (*it).second;
      }
   }
   return res;
}

/*----------------------------------------------------------------------------*/

unsigned char AaMap::_getId(const char* aa3) const {
   unsigned char id = AaMap::UNK;
   {
      /* returns 0 if not found, which corresponds to the aa code '?' */
      std::map<std::string, unsigned char>::const_iterator it = m->aa3map.find(
            aa3);
      if (it != m->aa3map.end()) {
         id = (*it).second;
      }
   }
   return id;
}

/*----------------------------------------------------------------------------*/

char AaMap::_getResn1(const char* resName3) const {
   char res1 = '?';
   {
      std::map<std::string, char>::const_iterator it = m->aa3toaa3map.find(
            resName3);
      if (it != m->aa3toaa3map.end()) {
         res1 = (*it).second;
      }
   }
   return res1;
}

/*----------------------------------------------------------------------------*/

char AaMap::_idxToResn1(uint32_t idx) const {
   char res = m->aa1table[0];
   {
      if (idx < m->aa1table.size()) {
         res = m->aa1table[idx];
      }
   }
   return res;
}

/*----------------------------------------------------------------------------*/

const char* AaMap::_getResn3(unsigned char resId) const {
   const char* res = NULL;
   {
      assert(resId < m->aa3table.size());
      res = m->aa3table[resId].c_str();
   }
   return res;
}

/*----------------------------------------------------------------------------*/

bool AaMap::_isAa(const char* aa3str) const {
   return (_getId(aa3str) != AaMap::UNK);
}

/*----------------------------------------------------------------------------*/

unsigned char AaMap::getId(char resName1) {
   return getInstance()._getId(resName1);
}

/*----------------------------------------------------------------------------*/

unsigned char AaMap::getId(const char* resName3) {
   return getInstance()._getId(resName3);
}

/*----------------------------------------------------------------------------*/

bool AaMap::isAa(const char* aa3str) {
   return getInstance()._isAa(aa3str);
}

/*----------------------------------------------------------------------------*/

char AaMap::getResn1(const char* resName3) {
   return getInstance()._getResn1(resName3);
}

/*----------------------------------------------------------------------------*/

char AaMap::idxToResn1(uint32_t idx) {
   return getInstance()._idxToResn1(idx);
}

/*----------------------------------------------------------------------------*/

const char* AaMap::getResn3(unsigned char resId) {
   return getInstance()._getResn3(resId);
}

/*----------------------------------------------------------------------------*/

AaMap::AaMap() {
   m = new AaMapPriv();
   init();
}

/*----------------------------------------------------------------------------*/

AaMap::~AaMap() {
   delete m;
}

/*----------------------------------------------------------------------------*/

AaMap& AaMap::getInstance() {
   static std::auto_ptr<AaMap> aamap;

   #pragma omp critical (aamap_init)
   if (aamap.get() == NULL) {
      aamap.reset(new AaMap());
   }
   return *aamap.get();
}
