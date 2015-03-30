#include "chainmatchexpansion.h"
#include <set>
#include <algorithm>

/*----------------------------------------------------------------------------*/

ChainMatchExpansion::ChainMatchExpansion() {
}

/*----------------------------------------------------------------------------*/

void ChainMatchExpansion::expandB2(ChainConnectionInfo &ci, SeqIdInfo &seqInfo,
         const Protein &pBC, const Protein &pU1, Chain *cB2,
         const ChainMatch &cm, DebugOut &dbg, const ConfigCompUnbound cfg) {

   const std::set<Chain*> sCmU1(cm.U1().begin(), cm.U1().end());

   m_B1.assign(cm.B1().begin(), cm.B1().end());

   /* which chains are connected to B2? */
   dbg("-- try to add chains from BX(%s) to B2 if connected", pBC.name().c_str());
   dbg("-- fails if similar to any of U1(%s)", pU1.name().c_str());
   /* B1: add matched B1 chains */
   std::set<Chain*> sB1(cm.B1().begin(), cm.B1().end());
   std::set<Chain*> sBX(cm.BX().begin(), cm.BX().end());
   /* BX: remove initial B2 chain */
   sBX.erase(cB2);
   std::set<Chain*> sB2;
   std::set<Chain*> sBS;
   /* B2: add initial chain */
   sB2.insert(cB2);

   dbg("try all connected chains (BX -> B2):");
   /* try to transfer chains of BX to B2 (require connectivity with B2) */
   while (sBX.empty() == false) {
      std::set<Chain*>::iterator sBXit;
      std::set<Chain*>::iterator sB2it;
      Chain* cx = NULL;
      /* try all of BX */
      for (sBXit = sBX.begin(); cx == NULL && sBXit != sBX.end(); ++sBXit) {
         /* search for connection to any of B2 */
         for (sB2it = sB2.begin(); sB2it != sB2.end(); ++sB2it) {
            /* has connection */
            const uint32_t intSize = ci.getInterfaceSize(*sB2it, *sBXit, *cfg.cmIntfDefinition);
            dbg("  %c -> %c : %3u intsize", (*sBXit)->name(),
                     (*sB2it)->name(), intSize);
            if (intSize >= cfg.minIntSize) {
               cx = (*sBXit);
               break;
            }
         }
      }
      if (cx != NULL) {
         /* similar to any of U1? */
         if (hasNoSimilar(cx, sCmU1, seqInfo, cfg.u1b2SeqIdMax, dbg) == true) {
            dbg("    added to B2");
            sB2.insert(cx);
         } else {
            dbg("    not added to B2! has similar in U1");
            sBS.insert(cx);
         }
         sBX.erase(cx);
      } else {
         break;
      }
   }
   m_BX.assign(sBX.begin(), sBX.end());
   m_B2.assign(sB2.begin(), sB2.end());
   m_BS.assign(sBS.begin(), sBS.end());
}

/*----------------------------------------------------------------------------*/

void ChainMatchExpansion::expandU1(ChainConnectionInfo &ci, SeqIdInfo &seqInfo,
         const Protein &pBC, const Protein &pU1, Chain *cB2,
         const ChainMatch &cm, DebugOut &dbg, const ConfigCompUnbound cfg) {
   /* get chains */
   const std::set<Chain*> sCmU1(cm.U1().begin(), cm.U1().end());

   m_B1.assign(cm.B1().begin(), cm.B1().end());

   /* which chains are connected to U1? */

   /* UX: add unmatched pU1 chains */
   std::set<Chain*> sUX(cm.UX().begin(), cm.UX().end());
   std::set<Chain*> sB2(m_B2.begin(), m_B2.end());
   std::set<Chain*> sU1(cm.U1().begin(), cm.U1().end());
   std::set<Chain*> sUC;

   /* try to transfer UX to U1
    * is not allowed to have an interface with any of B2
    * is required to have an interface with any of U1 */
   while (sUX.empty() == false) {
      std::set<Chain*>::iterator sUXit;
      std::set<Chain*>::iterator sU1it;
      std::set<Chain*>::iterator sB2it;
      Chain* cx = NULL;
      bool conflict = false;

      /* connected to U1? */
      for (sUXit = sUX.begin(); cx == NULL && sUXit != sUX.end(); ++sUXit) {
         /* search for connection to any of U1 */
         for (sU1it = sU1.begin(); sU1it != sU1.end(); ++sU1it) {
            /* has connection */
            const uint32_t intSize = ci.getInterfaceSize(*sU1it, *sUXit, *cfg.cmIntfDefinition);
            dbg("  ux->u1 %c %c intsize %3u", (*sUXit)->name(),
                     (*sU1it)->name(), intSize);
            if (intSize >= cfg.minIntSize) {
               cx = (*sUXit);
               break;
            }
         }
      }

      /* search for connection to any of B2 */
      if (cx != NULL) {
         for (sB2it = sB2.begin(); sB2it != sB2.end(); ++sB2it) {
            /* has connection/overlaps with B2? */
            const uint32_t intSize = ci.getInterfaceSize(*sB2it, cx, *cfg.cmIntfDefinition);
            dbg("  ux->b2 %c %c intsize %3u", cx->name(), (*sB2it)->name(),
                     intSize);
            if (intSize > cfg.minIntSize) {
               dbg("    conflict!");
               conflict = true;
               break;
            }
         }
      }
      if (cx != NULL) {
         if (conflict == true) {
            dbg("    add UC");
            sUC.insert(cx);
         } else {
            dbg("    add U1");
            sU1.insert(cx);
         }
         sUX.erase(cx);
      } else {
         /* no new candidate found - let's get out of here */
         break;
      }
   }
   m_U1.assign(cm.U1().begin(), cm.U1().end());
   /* add only new chains - to keep initial chains ordered */
   std::set_difference(sU1.begin(), sU1.end(), sCmU1.begin(), sCmU1.end(),
            std::back_inserter(m_U1));
   m_UX.assign(sUX.begin(), sUX.end());
   m_UC.assign(sUC.begin(), sUC.end());
}

/*----------------------------------------------------------------------------*/

void ChainMatchExpansion::fixUCligand(DebugOut &dbg, const ConfigCompUnbound cfg) {
   /* check number of chains */
   uint32_t numConflictingChains = m_UC.size();
   dbg("number of conflicting chains  : %u", numConflictingChains);
   if (numConflictingChains != 1) {
      return;
   }
   Chain* cc = m_UC.front();

   /* check number of residues */
   uint32_t numConflictingResis = cc->resis().size();
   dbg("number of residues of chain %c : %u/%u", cc->name(), numConflictingResis, cfg.cmxUCfixMaxResis);
   if (numConflictingResis > cfg.cmxUCfixMaxResis) {
      return;
   }

   /* check number of S2 bonds */
   uint32_t numS2Bonds = Chain::getS2Bonds(m_UC.front()->resis(), 2.3f).size();
   dbg("number of S2bonds of chain %c : %u/%u", cc->name(), numS2Bonds, cfg.cmxUCfixMaxS2bonds);
   if (numS2Bonds > cfg.cmxUCfixMaxS2bonds) {
      return;
   }

   /* transfer chain to solved */
   m_US.push_back(m_UC.front());
   m_UC.pop_back();
}

/*----------------------------------------------------------------------------*/

void ChainMatchExpansion::expandU1oligo(ChainConnectionInfo &ci,
         SeqIdInfo &seqInfo, const Protein &pBC, const Protein &pU1, Chain *cB2,
         const ChainMatch &cm, DebugOut &dbg, const ConfigCompUnbound cfg) {

   const std::set<Chain*> sCmU1(cm.U1().begin(), cm.U1().end());

   m_B1.assign(cm.B1().begin(), cm.B1().end());

   /* UX: add unmatched pU1 chains */
   std::set<Chain*> sUX(cm.UX().begin(), cm.UX().end());
   std::set<Chain*> sB2(m_B2.begin(), m_B2.end());
   std::set<Chain*> sU1(cm.U1().begin(), cm.U1().end());
   std::set<Chain*> sUC;
   std::set<Chain*> sUS;

   struct {
      uint32_t intSize; /* interface to U1 */
      Chain *c;
      bool conflict;
   } candRes;

   /* try to transfer UX to U1
    * is not allowed to have an interface with any of B2
    * is required to have an interface with any of U1
    * is not allowed to be seq. similar to any if U1 */

   while (sUX.empty() == false) {
      std::set<Chain*>::iterator sUXit;
      std::set<Chain*>::iterator sU1it;
      std::set<Chain*>::iterator sB2it;

      candRes.conflict = false;
      candRes.c = NULL;
      candRes.intSize = 0;

      /* search in UX for greatest connection to any of U1 */
      for (sUXit = sUX.begin(); sUXit != sUX.end(); ++sUXit) {
         for (sU1it = sU1.begin(); sU1it != sU1.end(); ++sU1it) {
            /* has connection */
            const uint32_t intSize = ci.getInterfaceSize(*sU1it, *sUXit, *cfg.cmIntfDefinition);
            dbg("  ux->u1 %c %c intsize %3u", (*sUXit)->name(),
                     (*sU1it)->name(), intSize);
            if (intSize >= cfg.minIntSize && intSize > candRes.intSize) {
               candRes.c = (*sUXit);
               candRes.intSize = intSize;
            }
         }
      }

      /* search for connection to any of B2 */
      if (candRes.c != NULL) {
         for (sB2it = sB2.begin(); sB2it != sB2.end(); ++sB2it) {
            /* has connection/overlaps with B2? */
            const uint32_t intSize = ci.getInterfaceSize(*sB2it, candRes.c, *cfg.cmIntfDefinition);
            dbg("  ux->b2 %c %c intsize %3u", candRes.c->name(),
                     (*sB2it)->name(), intSize);
            if (intSize > cfg.minIntSize) {
               dbg("    conflict!");
               candRes.conflict = true;
               break;
            }
         }
      }
      if (candRes.c != NULL) {
         /* only thread as conflict if not similar to any of U1 */
         if (hasNoSimilar(candRes.c, sCmU1, seqInfo, cfg.u1HomodimerSeqidMin, dbg) == true) {


            /* conflict ? */
            if (candRes.conflict) {
               dbg("    add UC");
               sUC.insert(candRes.c);
            } else {
               dbg("    add U1");
               sU1.insert(candRes.c);
            }
         } else {
            dbg("    add US");
            sUS.insert(candRes.c);
         }
         sUX.erase(candRes.c);
      } else {
         /* no new candidate found - let's get out of here */
         break;
      }
   }
   m_U1.assign(cm.U1().begin(), cm.U1().end());
   /* add only new chains - to keep initial chains ordered */
   std::set_difference(sU1.begin(), sU1.end(), sCmU1.begin(), sCmU1.end(),
            std::back_inserter(m_U1));
   m_UX.assign(sUX.begin(), sUX.end());
   m_UC.assign(sUC.begin(), sUC.end());
   m_US.assign(sUS.begin(), sUS.end());

}
