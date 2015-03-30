#ifndef CHAINOVERLAPINFO_H_
#define CHAINOVERLAPINFO_H_

#include "../libxtaldata/chain.h"
#include <stdint.h>
#include <map>

/*----------------------------------------------------------------------------*/

/* compute overlap from chainA to ChainB by counting the number of clashing
 * ca atoms (dist < minCaDist). finally divide by numCaAtoms(chainA) */

class ChainOverlapInfo {
public:

   ChainOverlapInfo(float minCaDist);

   float getChainOverlap(const Chain* chainB1, const Chain* chainB2,
         const Chain* chainU1);

private:
   std::map<std::pair<const Chain*, const Chain*>, float> m_overlap;

   float m_minCaDist;

private:
};

#endif /* CHAINOVERLAPINFO_H_ */
