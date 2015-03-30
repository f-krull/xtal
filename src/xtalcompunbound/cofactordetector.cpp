#include "cofactordetector.h"
#include "../libxtalcommon/chaininterfacetable.h"
#include <stdlib.h>
#include <stdio.h>
/*----------------------------------------------------------------------------*/

bool CofactorDetector::checkDna(const Protein &complex, const Chains &cB1, const Chains &cB2, const IntfDef &id, float minDnaDist) {

   bool hasDna = false;

   ChainInterfaceTable iftab(id);

   /* analyze all interfaces between the selected chains */
   for (uint32_t i = 0; i < cB1.size(); i++) {
      for (uint32_t j = 0; j < cB2.size(); ++j) {
         const ChainInterface *intf = iftab.get(cB1[i], cB2[j]);
         /* this will happen if chain c is not one of the two pot. intf. chains */
         if (intf == NULL) {
            continue;
         }
         /* which non-residues are in range of chain1s interface residues? */
         for (uint32_t k = 0; k < intf->int1.size(); k++) {
            if (intf->int1[k] == false) {
               /* not an interface residue -> continue */
               continue;
            }
            /* is any non-residue in range */
            for (uint32_t l = 0; l < complex.nonResis().size(); l++) {
               if (Residue::hasContact_fast(intf->c1->resis()[k],
                     complex.nonResis()[l], minDnaDist) == true) {
                  hasDna = hasDna || complex.nonResis()[l]->isNucleotide();
               }
            }
         }
         /* which non-residues are in range of chain1s interface residues? */
         for (uint32_t k = 0; k < intf->int2.size(); k++) {
            if (intf->int2[k] == false) {
               /* not an interface residue -> continue */
               continue;
            }
            /* is any non-residue in range */
            for (uint32_t l = 0; l < complex.nonResis().size(); l++) {
               if (Residue::hasContact(intf->c2->resis()[k],
                     complex.nonResis()[l], minDnaDist) == true) {
                  hasDna = hasDna || complex.nonResis()[l]->isNucleotide();
               }
            }
         }
      }
   }
   return (hasDna == false);
}

/*----------------------------------------------------------------------------*/

void setNeighbors(CofactorDetector::CofNeighborhood &cand,
         const Residues &interface, float minCofDist) {
   cand.neighbors.clear();

   for (uint32_t j = 0; j < interface.size(); j++) {
      if (Residue::hasContact(cand.cof, interface[j], minCofDist) == true) {
         /* add interface residue as neighbor */
         cand.neighbors.insert(interface[j]);
      }
   }
}


void CofactorDetector::check(const Protein &complex, const Residues &interface, float minCofDist, uint32_t minCofNeighbors, float cofNeighborDist, std::vector<std::string> &ignoreList, DebugOut &dbg) {
   m_hasDna = false;
   m_cofactors.clear();
   /* look at every non residue */
   for (uint32_t i = 0; i < complex.nonResis().size(); i++) {
      CofNeighborhood cand;
      cand.cof = complex.nonResis()[i];



      /* skip if HETATM has backbone */
      if (cand.cof->hasCompleteBackbone() == true) {
         continue;
      }


      bool ignore = false;
      for (uint32_t i = 0; i < ignoreList.size(); ++i) {
         /* skip leading blanks */
         const char* resnTrim = cand.cof->atoms().front()->resName();
         while (resnTrim[0] == ' ') {
            resnTrim++;
         }
         /* do not care about water, etc. */
         if (ignoreList[i] == resnTrim) {
            ignore = true;
            continue;
         }
      }
      if (ignore == true) {
         continue;
      }

      setNeighbors(cand, interface, minCofDist);

      /* candidate is nucleotide and within interface */
      if (cand.neighbors.size() > 0 && cand.cof->isNucleotide() == true) {
         m_hasDna = true;
      }

      if (std::string("ZNH") == cand.cof->atoms().front()->resName()) {
        dbg("ZHN num neigh %u", cand.neighbors.size());
     }

      /* see if it has at least x neighboring interface residues */
      if (cand.neighbors.size() >= minCofNeighbors) {
         /* expand neighbors at a larger cutoff for later */
         setNeighbors(cand, interface, cofNeighborDist);
         m_cofactors.push_back(cand);
      }
   }

   for (uint32_t i = 0; i < m_cofactors.size(); i++) {
      dbg("cofdetection: %s %s%c  '%s' numNeigh.=%u", complex.name().c_str(), m_cofactors[i].cof->atoms().front()->resSeq(), m_cofactors[i].cof->atoms().front()->chainId(),  m_cofactors[i].cof->atoms().front()->resName(),  m_cofactors[i].neighbors.size());
   }
   return;
}

/*----------------------------------------------------------------------------*/

CofactorDetector::CofactorDetector(const Protein &complex,
      const Residues &interface, float minCofDist, uint32_t minCofNeighbors, float cofNeighborDist, std::vector<std::string> &ignoreList, DebugOut &dbg) {
   m_hasDna = false;
   check(complex, interface, minCofDist, minCofNeighbors, cofNeighborDist, ignoreList, dbg);
}

/*----------------------------------------------------------------------------*/

const std::vector<CofactorDetector::CofNeighborhood> & CofactorDetector::cofactors() const {
   return m_cofactors;
}

/*----------------------------------------------------------------------------*/

bool CofactorDetector::hasDna() const {
   return m_hasDna;
}
