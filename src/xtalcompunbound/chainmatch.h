

#ifndef CHAINMATCH_H_
#define CHAINMATCH_H_

#include "../libxtaldata/protein.h"
#include "chainconnectioninfo.h"
#include "seqidinfo.h"
#include "configcompunbound.h"
#include "debugoutput.h"
#include <stdint.h>


class ChainMatch {
public:
   struct ChainPair {
      Chain *b; /* bound */
      Chain *u; /* unbound */
      bool operator<(const ChainPair &e) const {return b->name() < e.b->name();}
   };
   typedef std::vector<ChainPair> Matching;

   ChainMatch(ChainConnectionInfo &ci, SeqIdInfo &si, const Protein &pBC,
            Chain *cB1, Chain *cB2, const Protein &pU1, Chain *cU1, DebugOut &dbg);

   /* results */
   bool good() const;
   const Matching& getMatching() const;

   const Chains & B1() const; /* matched pBC chains */
   const Chains & BX() const; /* unmatched pBC chains */
   const Chains & U1() const; /* matched pU1 chains */
   const Chains & UX() const; /* unmatched pU1 chains */

private:
   Matching m_matching;

   Chains m_B1;
   Chains m_U1;
   Chains m_BX;
   Chains m_UX;

   bool m_good;
};

#endif /* CHAINMATCH_H_ */
