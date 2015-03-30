#include "atom.h"
#include "aamap.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <memory.h>
#include <string.h>

/*----------------------------------------------------------------------------*/

const char *Atom::RES_ALA = "ALA";
const char *Atom::RES_CYS = "CYS";
const char *Atom::RES_ASP = "ASP";
const char *Atom::RES_GLU = "GLU";
const char *Atom::RES_PHE = "PHE";
const char *Atom::RES_GLY = "GLY";
const char *Atom::RES_HIS = "HIS";
const char *Atom::RES_ILE = "ILE";
const char *Atom::RES_LYS = "LYS";
const char *Atom::RES_LEU = "LEU";
const char *Atom::RES_MET = "MET";
const char *Atom::RES_ASN = "ASN";
const char *Atom::RES_PRO = "PRO";
const char *Atom::RES_GLN = "GLN";
const char *Atom::RES_ARG = "ARG";
const char *Atom::RES_SER = "SER";
const char *Atom::RES_THR = "THR";
const char *Atom::RES_VAL = "VAL";
const char *Atom::RES_TRP = "TRP";
const char *Atom::RES_TYR = "TYR";
const char *Atom::RES_MSE = "MSE";
const char *Atom::RES_PYL = "PYL";
const char *Atom::RES_SEC = "SEC";
const char *Atom::RES_HSD = "HSD";
const char *Atom::RES_HSE = "HSE";
const char *Atom::RES_HSP = "HSP";
const char *Atom::RES_HOH = "HOH";

const char *Atom::RES_DA = " DA";
const char *Atom::RES_DC = " DC";
const char *Atom::RES_DG = " DG";
const char *Atom::RES_DT = " DT";

const char *Atom::RES_A = "  A";
const char *Atom::RES_C = "  C";
const char *Atom::RES_G = "  G";
const char *Atom::RES_T = "  T";

const char *Atom::ATM_CA  = " CA ";
const char *Atom::ATM_C   = " C  ";
const char *Atom::ATM_N   = " N  ";
const char *Atom::ATM_O   = " O  ";
const char *Atom::ATM_CB  = " CB ";
const char *Atom::ATM_P   = " P  ";
const char *Atom::ATM_O3p = " O3'";
const char *Atom::ATM_Sg  = " SG ";

const char Atom::ELMNT_C = 'C';
const char Atom::ELMNT_O = 'O';
const char Atom::ELMNT_H = 'H';
const char Atom::ELMNT_N = 'N';
const char Atom::ELMNT_S = 'S';
const char Atom::ELMNT_P = 'P';

/*----------------------------------------------------------------------------*/

#define READ_FIELD_STR(dst,src,offset,len) \
   if (offset + sizeof(dst)-1 < len) { \
      memcpy(dst, &src[offset], sizeof(dst)-1); \
      dst[sizeof(dst)-1] = '\0'; \
   } else { \
      memset(dst, 0, sizeof(dst)); \
   }

#define READ_FIELD_CHAR(dst,src,offset,len) \
   if (offset < len) { \
      dst = src[offset]; \
   } else {  \
      dst = '\0'; \
   }

#define READ_FIELD_FLOAT(dst,src,offset,len) \
   if (offset + 8 < len) { \
      dst = atof(&src[offset]); \
   } else { \
      dst = 0; \
   }

/*----------------------------------------------------------------------------*/

bool Atom::parse(const char* str) {
   size_t len = strlen(str);

   if (len < 54) {
      return false;
   }
   if (strncmp(str, "ATOM  ", 6) != 0 && strncmp(str, "HETATM", 6) != 0) {
      return false;
   }
   READ_FIELD_STR(m_rec, str, 0, len);
   READ_FIELD_STR(m_serial, str, 6, len);
   READ_FIELD_STR(m_name, str, 12, len);
   READ_FIELD_CHAR(m_altLoc, str, 16, len);
   READ_FIELD_STR(m_resName, str, 17, len);
   READ_FIELD_CHAR(m_chainId, str, 21, len);
   READ_FIELD_STR(m_resSeq, str, 22, len);
   READ_FIELD_CHAR(m_iCode, str, 26, len);
   READ_FIELD_FLOAT(m_v.f[0], str, 30, len);
   READ_FIELD_FLOAT(m_v.f[1], str, 38, len);
   READ_FIELD_FLOAT(m_v.f[2], str, 46, len);
   READ_FIELD_FLOAT(m_occupancy, str, 54, len);
   READ_FIELD_FLOAT(m_tempFactor, str, 60, len);
   READ_FIELD_STR(m_element, str, 76, len);
   READ_FIELD_STR(m_charge, str, 78, len);
   return true;
}

/*----------------------------------------------------------------------------*/

Atom* Atom::getParsed(const char *str) {
   Atom *a;

   a = new Atom();
   a->parse(str);
   return a;
}

/*----------------------------------------------------------------------------*/

Atom & Atom::operator=(const Vector &v2) {
   memcpy(this->m_v.f, v2.v(), sizeof(m_v.f));
   return *this;
}

/*----------------------------------------------------------------------------*/

void Atom::print() const {
   printf("%s", m_rec);
   printf("%s", m_serial);
   printf(".");
   printf("%s", m_name);
   printf("%c", m_altLoc);
   printf("%s", m_resName);
   printf(".");
   printf("%c", m_chainId);
   printf("%s", m_resSeq);
   printf("%c", m_iCode);
   printf("...");
   printf("%8.3f", m_v.f[0]);
   printf("%8.3f", m_v.f[1]);
   printf("%8.3f", m_v.f[2]);
   printf("%6.2f", m_occupancy);
   printf("%6.2f", m_tempFactor);
   printf("..........");
   printf("%s", m_element);
   printf("%s", m_charge);
   printf(".\n");
}

/*----------------------------------------------------------------------------*/

bool Atom::isCa() const {
   return strncmp(m_name, ATM_CA, 4) == 0;
}

/*----------------------------------------------------------------------------*/

bool Atom::isBackBone() const {
   bool res = false;

   res = res || strncmp(m_name, ATM_CA, 4) == 0;
   res = res || strncmp(m_name, ATM_C, 4) == 0;
   res = res || strncmp(m_name, ATM_N, 4) == 0;
   res = res || strncmp(m_name, ATM_O, 4) == 0;
   return res;
}

/*----------------------------------------------------------------------------*/

bool Atom::isElmnt(char c) const {
   return m_name[1] == c;
}

/*----------------------------------------------------------------------------*/

bool Atom::isResn(const char *resn) const {
   assert(strlen(resn) == 3);
   return strncmp(m_resName, resn, 3) == 0;
}

/*----------------------------------------------------------------------------*/

bool Atom::isAminoacid() {
   return AaMap::isAa(m_resName);
}

/*----------------------------------------------------------------------------*/

bool Atom::isNuclotide() const {
    bool res = false;

    res = res || strncmp(m_resName, Atom::RES_DA, 4) == 0;
    res = res || strncmp(m_resName, Atom::RES_DC, 4) == 0;
    res = res || strncmp(m_resName, Atom::RES_DG, 4) == 0;
    res = res || strncmp(m_resName, Atom::RES_DT, 4) == 0;
    res = res || strncmp(m_resName, Atom::RES_A, 4) == 0;
    res = res || strncmp(m_resName, Atom::RES_C, 4) == 0;
    res = res || strncmp(m_resName, Atom::RES_G, 4) == 0;
    res = res || strncmp(m_resName, Atom::RES_T, 4) == 0;
    return res;
}

/*----------------------------------------------------------------------------*/

bool Atom::isHetatm() const {
   return m_rec[0] == 'H';
}

/*----------------------------------------------------------------------------*/

Vector Atom::centre(const Atoms & atoms)
{
    Vector centre;
    for(unsigned int i = 0;i < atoms.size();i++){
        centre = *atoms[i] + centre;
    }
    return centre /= atoms.size();
}

/*----------------------------------------------------------------------------*/

bool Atom::isName(const char *name) const {
   assert(strlen(name) == 4);
   return strncmp(m_name, name, 4) == 0;
}
