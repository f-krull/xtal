#ifndef CHAINMATCHEXPANSION_H_
#define CHAINMATCHEXPANSION_H_

/*----------------------------------------------------------------------------*/

#include "chainmatch.h"
#include "ichainmatchexp.h"
#include "debugoutput.h"

/*----------------------------------------------------------------------------*/

class ChainMatchExpansion: public IChainMatchExp {
public:
   ChainMatchExpansion();

   /* can be done before we computed any interface */
   void expandB2(ChainConnectionInfo &ci, SeqIdInfo &seqInfo,
            const Protein &pBC, const Protein &pU1, Chain *cB2,
            const ChainMatch &cm, DebugOut &dbg, const ConfigCompUnbound cfg =
                     ConfigCompUnbound());

   /* pU needs to be superimposed already */
   void expandU1(ChainConnectionInfo &ci, SeqIdInfo &seqInfo,
            const Protein &pBC, const Protein &pU1, Chain *cB2,
            const ChainMatch &cm, DebugOut &dbg, const ConfigCompUnbound cfg =
                     ConfigCompUnbound());

   /* if UC consists of one single - allow it if it has no S2-bonds */
   void fixUCligand(DebugOut &dbg, const ConfigCompUnbound cfg = ConfigCompUnbound());


   /* fall-back option if we observed any UC */
   void expandU1oligo(ChainConnectionInfo &ci, SeqIdInfo &seqInfo,
            const Protein &pBC, const Protein &pU1, Chain *cB2,
            const ChainMatch &cm, DebugOut &dbg, const ConfigCompUnbound cfg =
                     ConfigCompUnbound());


private:
};

#endif /* CHAINMATCHEXPANSION_H_ */
