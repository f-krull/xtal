

#ifndef COFACTORDETECTOR_H_
#define COFACTORDETECTOR_H_

#include "../libxtaldata/protein.h"
#include "../libxtalcommon/interface.h"
#include "debugoutput.h"
#include <set>
#include <stdint.h>


/*----------------------------------------------------------------------------*/
/* detects cofactors and provides assigned neighboring residues */
class CofactorDetector {
public:

   struct CofNeighborhood {
      Residue *cof;
      std::set<Residue*> neighbors;
   };

   CofactorDetector(const Protein &complex, const Residues &interface,
         float minCofDist, uint32_t minCofNeighbors, float cofNeighborDist, std::vector<std::string> &ignoreList, DebugOut &dbg);

   const std::vector<CofNeighborhood> &cofactors() const;
   bool hasDna() const;

   static bool checkDna(const Protein &complex, const Chains &cB1,
         const Chains &cB2, const IntfDef &id, float minDnaDist);

private:
   void check(const Protein &complex, const Residues &interface,
         float minCofDist, uint32_t minCofNeighbors, float cofNeighborDist, std::vector<std::string> &ignoreList, DebugOut &dbg);

   std::vector<CofNeighborhood> m_cofactors;
   bool m_hasDna;
};

#endif /* COFACTORDETECTOR_H_ */
