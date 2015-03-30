#ifndef IATOMS_H_
#define IATOMS_H_

#include "atom.h"
#include <vector>

/*----------------------------------------------------------------------------*/

class PdbInfo {
public:
   char title[41];
   char date[10];
   char id[5];
   char resolution[6];
   char expdata[61];

   PdbInfo() {
      title[0] = '\0';
      date[0] = '\0';
      id[0] = '\0';
      resolution[0] = '\0';
      expdata[0] = '\0';
   }
};

/*----------------------------------------------------------------------------*/

class IAtoms {
public:

   virtual ~IAtoms();

   const virtual Atoms & atoms() const;

   virtual Atoms & atoms();

   virtual size_t size() const;

   virtual bool readPdb(const char *filename, PdbInfo* pdbinf = NULL);

   virtual Atom* getByName(const char *name) const;

   Vector centre() const;

   void print();

protected:
   IAtoms();

   template<typename T> void set(T begin, T end);

   template<typename T> void insert(T begin, T end);

   void add(Atom *a);

   void destroy();

   virtual void update();

protected:
   Atoms m_atoms;
};



/*----------------------------------------------------------------------------*/

template<typename T> inline void IAtoms::set(T begin, T end)
{
    m_atoms.assign(begin, end);
    update();
}

#endif /* IATOMS_H_ */
