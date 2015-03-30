#include "chaininterfaceinfo.h"
#include "../libxtalcommon/intfseqcomparer.h"

/*----------------------------------------------------------------------------*/


float ChainInterfaceInfo::getIntSim(const Chain *b1, const Chain* b2,
      const Chain *u1, const Chain *u2) {
   const std::pair<ChainInt, ChainInt> key = std::make_pair(
         std::make_pair(b1, b2), std::make_pair(u1, u2));
   std::map<std::pair<ChainInt, ChainInt>, float>::const_iterator it =
         m_intSim.find(key);

   /* compute only if not found */
   if (it != m_intSim.end()) {
      return it->second;
   }
   /* prepare input */
   std::vector<const Chain*> cb1;
   cb1.push_back(b1);
   std::vector<const Chain*> cb2;
   cb2.push_back(b2);
   std::vector<const Chain*> cu1;
   cu1.push_back(u1);
   std::vector<const Chain*> cu2;
   cu2.push_back(u2);

   /* get sequence with highlighted resis */
   std::string ib1 = m_cIntfTab.getSeqInt(b1, cb2);
   std::string ib2 = m_cIntfTab.getSeqInt(b2, cb1);
   std::string iu1 = m_cIntfTab.getSeqInt(u1, cu2);
   std::string iu2 = m_cIntfTab.getSeqInt(u2, cu1);

   /* align both interfaces sides */
   IntfSeqComparer isCmp;
   isCmp.alignMin(ib1, iu1);
   const float score1 = isCmp.getIntScore();
   isCmp.alignMin(ib2, iu2);
   const float score2 = isCmp.getIntScore();
   /* combine */
   m_intSim[key] = IntfSeqComparer::combineIntSim(score1, score2);
   return m_intSim[key];
}

