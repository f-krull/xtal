

#ifndef CHAINCONNECTIONINFO_H_
#define CHAINCONNECTIONINFO_H_

#include "../libxtaldata/chain.h"
#include "../libxtalcommon/interface.h"
#include <stdint.h>
#include <map>

/*----------------------------------------------------------------------------*/

class ChainConnectionInfo {
public:

   ChainConnectionInfo();

   uint32_t getInterfaceSize(const Chain*, const Chain*, const IntfDef &id);

private:
   std::map<std::pair<const Chain*, const Chain*>, float> m_intSize;
};


#endif /* CHAINCONNECTIONINFO_H_ */
