#ifndef SEQIDINFO_H_
#define SEQIDINFO_H_


#include "../libxtaldata/chain.h"
#include <stdint.h>
#include <map>

/*----------------------------------------------------------------------------*/

class SeqIdInfo {
public:
   struct IdInf {
     float min;
     float max;
     bool valid;
   };

   SeqIdInfo();

   const IdInf& idInf(const Chain* cB, const Chain* cU, bool audoCompute = true);

private:

   std::map<std::pair<const Chain*, const Chain*>, IdInf> m_seqId;
   IdInf m_sm;
};

#endif /* CHAINSEQIDINFO_H_ */
