#ifndef IRESIDUES_H_
#define IRESIDUES_H_

#include "residue.h"

/*----------------------------------------------------------------------------*/

class IResidues {
public:
   virtual ~IResidues();

   const virtual Residues & resis() const;

   virtual Residues & resis();

   const virtual Residues & nonResis() const;

   virtual Residues & nonResis();

   virtual size_t size() const;

   virtual bool build(IAtoms &atoms);

   Vector centre() const;

   void calcBonds();
   void calcSses();

protected:
   IResidues();

   template<typename T> void set(T begin, T end) {
      m_resis.assign(begin, end);
      update();
   }

   void add(Residue*);

   void destroy();

   virtual void update();

protected:
   Residues m_resis;
   Residues m_nonRes;
};


#endif /* IRESIDUES_H_ */
