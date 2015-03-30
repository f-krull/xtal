#include "chainmatch.h"
#include "stdint.h"
#include "chaininterfaceinfo.h"
#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <typeinfo>
#include <set>

/*----------------------------------------------------------------------------*/

struct ChainSetItem {
   Chain* c;
   bool valid;
};

typedef std::vector<ChainSetItem> ChainSet;

struct MatchInfo {
   const ConfigCompUnbound &cfg;
   ChainConnectionInfo &ci;
   ChainInterfaceInfo &ii;
   SeqIdInfo &si_BU;
   const uint32_t numMaxSteps;
   const uint32_t numMaxMatched;
};

/*----------------------------------------------------------------------------*/

/* check if pair has similar sequence (node label) */
static inline bool check_sequence(const ChainMatch::ChainPair &p,
      const MatchInfo &mi) {
   bool ok = false;
   SeqIdInfo::IdInf idInf = mi.si_BU.idInf(p.b, p.u);
   /* check sequence identity */
   if (idInf.valid == true && idInf.max > mi.cfg.chainMatchMinSeqId) {
      ok = true;
   }
   return ok;
}

/*----------------------------------------------------------------------------*/

/* check if interfaces to every matched pair are similar (edge label) */
static bool check_connectivity(const ChainMatch::ChainPair &p,
      ChainMatch::Matching &vM, const MatchInfo &mi, const ConfigCompUnbound &cfg, DebugOut &dbg) {

   bool ret = false;

   /* we don't need to check connectivity, if matching list is empty */
   if (vM.size() == 0) {
      return true;
   }

   for (uint32_t i = 0; i < vM.size(); i++) {

      const int32_t intB = mi.ci.getInterfaceSize(p.b, vM[i].b, *cfg.cmIntfDefinition);
      const int32_t intU = mi.ci.getInterfaceSize(p.u, vM[i].u, *cfg.cmIntfDefinition);

      /* chain is bound to matching in ub? */
      if (intB == 0) {
         /*no... */
         continue;
      }

      /* at least one significant interface? */
      if (intB < (int32_t)mi.cfg.minIntSize && intU < (int32_t)mi.cfg.minIntSize) {
         /* no... */
         continue;
      }


      /* do interface differ significantly? */
      const float intSim = mi.ii.getIntSim(p.b, vM[i].b, p.u, vM[i].u);
      if (intSim > mi.cfg.cmMinInfSim) {
         /* we found the same connection in bound and unbound */
         //dbg("intsim (%c,%c)(%c,%c): %5.3f", p.b->name(), vM[i].b->name(), p.u->name(), vM[i].u->name(), );
         ret = true;
      } else {
         /* we found significant interface differences - bail out */
         return false;
      }
   }
   return ret;
}

/*----------------------------------------------------------------------------*/

static ChainSet chainsset_init(const Chains &c) {
   ChainSet cs;
   for (uint32_t i = 0; i < c.size(); i++) {
      cs.push_back(ChainSetItem());
      cs.back().c = c[i];
      cs.back().valid = true;
   }
   return cs;
}

/*----------------------------------------------------------------------------*/

static bool chainsset_remove(ChainSet &cs, const Chain *c) {
   for (uint32_t i = 0; i < cs.size(); i++) {
      if (cs[i].c == c) {
         cs[i].valid = false;
      }
   }
   return false;
}

/*----------------------------------------------------------------------------*/

/**
 * try to match all chains of unbound (pU1) to chains of bound (pB1)
 *
 * start with mapping M (cB1,cU1); blacklist B (cB2,x), candidates C (cBX,cUX)
 *
 * function map (M, C):
 * if |C| > 0 (candidates left):
 *    for all candidates of C (cBX,...) connected to mapping M (cBY,..)?
 *       if can be mapped to an unbound? (cBX, cUX)?
 *          (by sequence)
 *          (by similar connectivity to all previous mappings)
 *          M = M + (cBX, cUX)
 *          C = C - (cBX, cUX)
 *          map(M, C)
 * else (no candidates left)
 *    return M
 */
void match_rec(ChainSet &sCB1, ChainSet &sCU1, const ChainMatch::ChainPair &pin,
      ChainMatch::Matching &vM, const MatchInfo &mi, ChainMatch::Matching &vR, uint32_t *numOp, const ConfigCompUnbound &cfg, DebugOut &dbg) {
   (*numOp)++;
   if ((*numOp) % 1000000 == 0) {
      dbg("numOp %u", *numOp);
   }

   //Log::dbg("%*s(%c,%c)", vM.size() * 2, "", pin.b->name(), pin.u->name());
   if (check_sequence(pin, mi) == false) {
      return;
   }
   //Log::dbg("%*sseq %s", vM.size() * 2 + 1, "", "Ok");
   if (check_connectivity(pin, vM, mi, cfg, dbg) == false) {
      return;
   }
   //Log::dbg("%*sint %s", vM.size() * 2 + 1, "", "Ok");

   vM.push_back(pin);

   /* explore new matches */
   ChainMatch::ChainPair pout;
   for (uint32_t i = 0; i < sCU1.size(); i++) {
      if (sCU1[i].valid == false) {
         continue;
      }
      pout.u = sCU1[i].c;
      sCU1[i].valid = false;
      for (uint32_t j = 0; j < sCB1.size(); j++) {

         if (sCB1[j].valid == false) {
            continue;
         }
         /* bail out at limit */
         if (*numOp > mi.numMaxSteps) {
            return;
         } else if (vR.size() >= mi.numMaxMatched) {
            return;
         }
         pout.b = sCB1[j].c;
         sCB1[j].valid = false;
         match_rec(sCB1, sCU1, pout, vM, mi, vR, numOp, cfg, dbg);
         sCB1[j].valid = true;
      }
      sCU1[i].valid = true;
   }
   /* found a new result? */
   if (vM.size() > vR.size()) {
      vR = vM;
   }

   vM.pop_back();
}

/*----------------------------------------------------------------------------*/
/* TODO move check to caller class */
static Chains shortChainsFilter(const Chains &chains, uint32_t minLen) {
   Chains res;
   for (uint32_t i = 0; i < chains.size(); i++) {
      uint32_t numStdResis = 0;
      for (uint32_t j = 0; j < chains[i]->resis().size(); j++) {
         if (chains[i]->resis()[j]->hasStandardCa()) {
            numStdResis++;
         }
      }
      if (numStdResis >= minLen) {
         res.push_back(chains[i]);
      }
   }
   return res;
}

/*----------------------------------------------------------------------------*/

#include "ubhelper.h"

ChainMatch::ChainMatch(ChainConnectionInfo &ci, SeqIdInfo &si, const Protein &pBC, Chain *cB1, Chain *cB2, const Protein &pU1,
      Chain *cU1, DebugOut &dbg) {
   m_good = false;
   /* good initial match ? - all chains are translated? */
   if (cU1 == NULL) {
      return;
   }
   if (cB1 == NULL) {
      return;
   }

   const ConfigCompUnbound cfg;

   /* get chains */
   Chains cCB1 = shortChainsFilter(pBC.chains(), cfg.minSeqLen);
   Chains cCU1 = shortChainsFilter(pU1.chains(), cfg.minSeqLen);
   dbg("chain sets (after filtering short chains):");
   dbg("  %lu BC chains (%s)", cCB1.size(), chainsToStr(cCB1).c_str());
   dbg("  %lu U1 chains (%s)", cCU1.size(), chainsToStr(cCU1).c_str());

   /* prebuild matching info */

   const uint32_t numMaxSteps = 100000000;
   const uint32_t numMaxMatchted = std::min(cCU1.size(), cCB1.size());
   ChainInterfaceInfo ii(*cfg.cmIntfDefinition);
   const MatchInfo mi = { cfg, ci, ii, si, numMaxSteps, numMaxMatchted };
   ChainSet sCB1 = chainsset_init(cCB1);
   ChainSet sCU1 = chainsset_init(cCU1);


   Matching vM;
   Matching &vR = m_matching;

   /* add initial mapping */
   ChainPair p;
   p.b = cB1;
   p.u = cU1;

   /* check for initial matching has to be done by aoller */
   /* ... */

   /* remove initial matching from candidate sets */
   chainsset_remove(sCB1, p.b);
   chainsset_remove(sCU1, p.u);

   /* match */
   uint32_t numOp = 0;
   match_rec(sCB1, sCU1, p, vM, mi, vR, &numOp, cfg, dbg);

   dbg("computed pairings:");
   for (uint32_t i = 0; i < vR.size(); i++) {
      dbg("    %u (%c,%c)", i, vR[i].b->name(), vR[i].u->name());
   }
   /* we found at least one match? */
   m_good = vR.size() > 0;

   /* matched chains */
   m_B1.reserve(m_matching.size());
   m_U1.reserve(m_matching.size());
   for (uint32_t i = 0; i < m_matching.size(); i++) {
      m_B1.push_back(m_matching[i].b);
      m_U1.push_back(m_matching[i].u);
   }


   std::sort(cCB1.begin(), cCB1.end());
   std::sort(cCU1.begin(), cCU1.end());
   std::sort(m_B1.begin(), m_B1.end());
   std::sort(m_U1.begin(), m_U1.end());
   /* put unmatched chains of pB to BX */
   std::set_difference(cCB1.begin(), cCB1.end(), m_B1.begin(), m_B1.end(),
            std::inserter(m_BX, m_BX.begin()));
   /* put unmatched chains of pU to UX */
   std::set_difference(cCU1.begin(), cCU1.end(), m_U1.begin(), m_U1.end(),
            std::inserter(m_UX, m_UX.begin()));

   /* make sure we have not matched cB2 (has to appear in m_BX) */
   m_good = m_good && std::find(m_BX.begin(), m_BX.end(), cB2) != m_BX.end();

   /* sort matchings and reassign B1/U1 (got sorted independently) */
   m_B1.clear();
   m_U1.clear();
   std::sort(m_matching.begin(), m_matching.end());
   for (uint32_t i = 0; i < m_matching.size(); i++) {
      m_B1.push_back(m_matching[i].b);
      m_U1.push_back(m_matching[i].u);
   }

   dbg("computed seq.IDs (%s(B1) -> %s(U1))", pBC.name().c_str(), pU1.name().c_str());
   for (uint32_t i = 0; i < cCB1.size(); i++) {
      for (uint32_t j = 0; j < cCU1.size(); j++) {
         const SeqIdInfo::IdInf &idi = si.idInf(cCB1[i], cCU1[j], true);
         if (idi.valid == true) {
            dbg("  chain %c -> %c : %.2f (seq.ID)", cCB1[i]->name(), cCU1[j]->name(), idi.max);
         }
      }
   }

   dbg("match stat %2lu/%2lu B1:%lu U1:%lu B2:%lu", cCB1.size(), cCU1.size(), B1().size(), U1().size(), BX().size());

   uint32_t ind = 2;
   dbg("chain match: ");
   dbg("%*s B1: %s", ind, "", chainsToStr(B1(), false).c_str());
   dbg("%*s BX: %s", ind, "", chainsToStr(BX()).c_str());
   dbg("%*s U1: %s", ind, "", chainsToStr(U1(), false).c_str());
   dbg("%*s UX: %s", ind, "", chainsToStr(UX()).c_str());
}

/*----------------------------------------------------------------------------*/

bool ChainMatch::good() const {
   return m_good;
}

/*----------------------------------------------------------------------------*/


const ChainMatch::Matching& ChainMatch::getMatching() const {
   return m_matching;
}


const Chains & ChainMatch::B1() const {
   return m_B1;
}

/*----------------------------------------------------------------------------*/

const Chains & ChainMatch::BX() const {
   return m_BX;
}

/*----------------------------------------------------------------------------*/

const Chains & ChainMatch::U1() const {
   return m_U1;
}

/*----------------------------------------------------------------------------*/

const Chains & ChainMatch::UX() const {
   return m_UX;
}
