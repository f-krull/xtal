#ifndef INTERFACE_H_
#define INTERFACE_H_

#include <stdint.h>

/*----------------------------------------------------------------------------*/

#include "../libxtaldata/residue.h"
#include <set>

typedef char cbool;


class IntfDef {
public:
   virtual ~IntfDef() {};
   virtual cbool operator()(const Residue *r, const Residue* o) const = 0;
};

/*----------------------------------------------------------------------------*/

class IntfDefAllAtom : public IntfDef {
public:
   IntfDefAllAtom(float cutoff);
   virtual ~IntfDefAllAtom() {};
   virtual cbool operator()(const Residue *r, const Residue* o) const;
private:
   float m_cutoff;
};

/*----------------------------------------------------------------------------*/

class IntfDefMendez2003 : public IntfDef {
public:
   virtual ~IntfDefMendez2003() {};
   virtual cbool operator()(const Residue *r, const Residue* o) const;
};

/*----------------------------------------------------------------------------*/


class Interface {
public:
   Interface(const std::vector<Residue*> &r, const std::vector<Residue*> &l,
            const IntfDef &id);

   ~Interface();

   std::vector<Residue*> getIntR() const;
   std::vector<Residue*> getIntL() const;

   float getInterfaceRmsd() const;

   static void getInterface(const std::vector<Residue*> &chain1,
            const std::vector<Residue*> &chain2, std::vector<cbool> &int1,
            std::vector<cbool> &int2, const IntfDef *id);

   static uint32_t getInterfaceSize(const Residues &r1, const Residues &r2,
            const IntfDef &id);

protected:
   std::set<Residue*> m_ir;
   std::set<Residue*> m_il;
   std::vector<Atom*> m_ica;
   std::vector<Atom*> m_ica_cpy;
};





#endif /* INTERFACE_H_ */
