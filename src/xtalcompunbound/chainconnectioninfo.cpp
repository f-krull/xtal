#include "chainconnectioninfo.h"
#include "../libxtalcommon/interface.h"
#include "../libxtalutil/log.h"

/*----------------------------------------------------------------------------*/

ChainConnectionInfo::ChainConnectionInfo() {
#if 0
   /* compute all pairwise interface sizes */
   for (uint32_t i = 0; i < chains.size(); i++) {
      for (uint32_t j = i + 1; j < chains.size(); j++) {
         uint32_t intSize = Interface::getInterfaceSize(chains[i]->resis(),
               chains[j]->resis());
         m_intSize[std::make_pair(chains[i], chains[j])] = intSize;
         m_intSize[std::make_pair(chains[j], chains[i])] = intSize;
      }
   }
   for (uint32_t i = 0; i < chains.size(); i++) {
      for (uint32_t j = i + 1; j < chains.size(); j++) {
         uint32_t intf = getInterfaceSize(chains[i], chains[j]);
         if (intf > 0) {
            Log::dbg("Interface %c:%c %u", chains[i]->name(), chains[j]->name(), getInterfaceSize(chains[i], chains[j]));
         }
      }
   }
#endif
}

/*----------------------------------------------------------------------------*/

uint32_t ChainConnectionInfo::getInterfaceSize(const Chain* chainA,
      const Chain* chainB, const IntfDef &id) {
   std::map<std::pair<const Chain*, const Chain*>, float>::const_iterator it;

   it = m_intSize.find(std::make_pair(chainA, chainB));
   /* not contained - return default */
   if (it == m_intSize.end()) {
      uint32_t intSize = Interface::getInterfaceSize(chainA->resis(),
            chainB->resis(), id);
      m_intSize[std::make_pair(chainB, chainA)] = intSize;
      m_intSize[std::make_pair(chainA, chainB)] = intSize;
      return intSize;
   }
   /* return interface size */
   return it->second;
}
