#include "chainoverlapinfo.h"

/*----------------------------------------------------------------------------*/

ChainOverlapInfo::ChainOverlapInfo(float minCaDist) {
   m_minCaDist = minCaDist;
}

/*----------------------------------------------------------------------------*/

static float computeChainOverlap(const Residues& chainA, const Residues& chainB, float minCaDist) {
   uint32_t numCa = 0;
   uint32_t numOverlap = 0;
   for (uint32_t i = 0; i < chainA.size(); i++) {
      if (chainA[i]->hasCa() == false) {
         continue;
      }
      numCa++;
      for (uint32_t j = 0; j < chainB.size(); j++) {
         if (chainB[j]->hasCa() == false) {
            continue;
         }
         if (chainB[j]->ca()->dist(*chainA[i]->ca()) < minCaDist) {
            /* overlap found - use next Ca of ChainA */
            numOverlap++;
            break;
         }
      }
   }

   return float(numOverlap) / numCa;
}

/*----------------------------------------------------------------------------*/


float ChainOverlapInfo::getChainOverlap(const Chain* chainB1, const Chain* chainB2,
      const Chain* chainU1) {
   std::map<std::pair<const Chain*, const Chain*>, float>::const_iterator it;

   it = m_overlap.find(std::make_pair(chainU1, chainB2));
   /* not contained - return default */
   if (it == m_overlap.end()) {
      float overlap = computeChainOverlap(chainB2->resis(), chainU1->resis(),
            m_minCaDist);
      m_overlap[std::make_pair(chainU1, chainB2)] = overlap;
      return overlap;
   }
   /* return interface size */
   return it->second;
}

