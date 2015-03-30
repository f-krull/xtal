#ifndef PROTEIN_H_
#define PROTEIN_H_

#include "../libxtaldata/mol.h"
#include "iresidues.h"
#include "ichains.h"

/*----------------------------------------------------------------------------*/

class Protein : public Mol, public IResidues, public IChains {
public:

   Protein();
   Protein(const Protein &p);
   virtual ~Protein();

   virtual bool readPdb(const char *filename);
   void info();
   Vector centre() const;

protected:

};

#endif /* PROTEIN_H_ */

