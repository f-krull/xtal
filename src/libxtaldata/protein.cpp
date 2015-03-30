#include "protein.h"
#include <stdio.h>
#include <stdint.h>
#include "../libxtalutil/log.h"

/*----------------------------------------------------------------------------*/

Protein::~Protein() {
   IResidues::destroy();
   IChains::destroy();
}

Protein::Protein()
{
}

/*----------------------------------------------------------------------------*/

Protein::Protein(const Protein & p) :
      Mol(p) {
   IResidues::build(*this);
   IChains::build(*this,*this);
}

/*----------------------------------------------------------------------------*/

bool Protein::readPdb(const char *filename) {
   bool ok = true;

   ok = ok && Mol::readPdb(filename);
   ok = ok && IResidues::build(*this);
   ok = ok && IChains::build(*this,*this);
   return ok;
}

/*----------------------------------------------------------------------------*/

void Protein::info() {
   Log::inf("Protein %s    chains:%lu resis:%lu atoms:%lu    nonresis:%lu",
         this->name().c_str(), chains().size(), resis().size(), atoms().size(),
         nonResis().size());
   for (uint32_t i = 0; i < chains().size(); i++) {
      Log::inf("  chain %c residues %lu %lu", chains()[i]->name(),
            chains()[i]->resis().size(), chains()[i]->atoms().size());
   }
}

/*----------------------------------------------------------------------------*/

Vector Protein::centre() const
{
    return IChains::centre();
}


