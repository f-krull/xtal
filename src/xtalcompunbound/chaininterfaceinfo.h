#ifndef CHAININTERFACEINFO_H_
#define CHAININTERFACEINFO_H_


#include "../libxtaldata/chain.h"
#include "../libxtalcommon/chaininterfacetable.h"
#include <map>

/*----------------------------------------------------------------------------*/

/* compare two interfaces - uses interface residue coverage */
class ChainInterfaceInfo {
public:

   typedef std::pair<const Chain*, const Chain*> ChainInt;

   ChainInterfaceInfo(const IntfDef &id) : m_cIntfTab(id) {
   }

   float getIntSim(const Chain *b1, const Chain* b2, const Chain *u1,
         const Chain *u2);

private:
   std::map<std::pair<ChainInt, ChainInt>, float> m_intSim;

   ChainInterfaceTable m_cIntfTab;
};


#endif /* CHAININTERFACEINFO_H_ */
