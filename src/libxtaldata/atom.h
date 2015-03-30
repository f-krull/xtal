
#ifndef ATOM_H_
#define ATOM_H_

#include "vector.h"
#include <vector>

/*----------------------------------------------------------------------------*/

class Atom;

/*----------------------------------------------------------------------------*/

typedef std::vector<Atom*> Atoms;

/*----------------------------------------------------------------------------*/

class Atom : public Vector {
public:

   void print() const;

   bool parse(const char* str);

   static Atom* getParsed(const char* str);

   const char* resName() const {
      return m_resName;
   }

   const char* resSeq() const {
      return m_resSeq;
   }

   const char iCode() const {
      return m_iCode;
   }

   const char* name() const {
      return m_name;
   }

   char chainId() const {
      return m_chainId;
   }

   float tempFactor() const {
      return m_tempFactor;
   }

   bool isCa() const;

   bool isBackBone() const;

   bool isName(const char *name) const;

   bool isResn(const char* resn) const;

   bool isElmnt(char type) const;

   bool isNuclotide() const;

   bool isAminoacid();

   bool isHetatm() const;

   Atom & operator=(const Vector &v2);

   /* residues */
   static const char* RES_ALA;
   static const char* RES_CYS;
   static const char* RES_ASP;
   static const char* RES_GLU;
   static const char* RES_PHE;
   static const char* RES_GLY;
   static const char* RES_HIS;
   static const char* RES_ILE;
   static const char* RES_LYS;
   static const char* RES_LEU;
   static const char* RES_MET;
   static const char* RES_ASN;
   static const char* RES_PRO;
   static const char* RES_GLN;
   static const char* RES_ARG;
   static const char* RES_SER;
   static const char* RES_THR;
   static const char* RES_VAL;
   static const char* RES_TRP;
   static const char* RES_TYR;
   static const char* RES_MSE;
   static const char* RES_PYL;
   static const char* RES_SEC;
   static const char* RES_HSD;
   static const char* RES_HSE;
   static const char* RES_HSP;

   /* DNA bases */
   static const char* RES_DA;
   static const char* RES_DC;
   static const char* RES_DG;
   static const char* RES_DT;
   /* RNA bases */
   static const char* RES_A;
   static const char* RES_C;
   static const char* RES_G;
   static const char* RES_T;

   /* non residues */
   static const char* RES_HOH;

   static const char* ATM_CA;
   static const char* ATM_C;
   static const char* ATM_O;
   static const char* ATM_N;
   static const char* ATM_CB;
   static const char* ATM_P;
   static const char* ATM_O3p;
   static const char* ATM_Sg;

   static const char ELMNT_C;
   static const char ELMNT_O;
   static const char ELMNT_H;
   static const char ELMNT_N;
   static const char ELMNT_S;
   static const char ELMNT_P;

   static Vector centre(const Atoms &atoms);

protected:
   char m_rec[7];
   char m_serial[6];
   char m_name[5];
   char m_altLoc;
   char m_resName[4];
   char m_chainId;
   char m_resSeq[5];
   char m_iCode;
   float m_occupancy;
   float m_tempFactor;
   char m_element[3];
   char m_charge[3];
};

#endif /* ATOM_H_ */


// 1 -  6        Record name   "ATOM  "
// 7 - 11        Integer       serial       Atom  serial number.
//13 - 16        Atom          name         Atom name.
//17             Character     altLoc       Alternate location indicator.
//18 - 20        Residue name  resName      Residue name.
//22             Character     chainID      Chain identifier.
//23 - 26        Integer       resSeq       Residue sequence number.
//27             AChar         iCode        Code for insertion of residues.
//31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
//39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
//47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
//55 - 60        Real(6.2)     occupancy    Occupancy.
//61 - 66        Real(6.2)     tempFactor   Temperature  factor.
//77 - 78        LString(2)    element      Element symbol, right-justified.
//79 - 80        LString(2)    charge       Charge  on the atom.

