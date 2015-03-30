#include "chainmatchexpansion.h"
#include "ubhelper.h"
#include "debugoutput.h"
#include <typeinfo>
#include <algorithm>

/*----------------------------------------------------------------------------*/

bool IChainMatchExp::hasNoSimilar(const Chain *cx, const std::set<Chain*> &sCmU1,
         SeqIdInfo &si, float maxSeqId, DebugOut &dbg) {

   std::set<Chain*>::iterator sCmU1it;

   for (sCmU1it = sCmU1.begin(); sCmU1it != sCmU1.end(); ++sCmU1it) {
      if (si.idInf(cx, *sCmU1it).max > maxSeqId) {
         return false;
      }
   }

   return true;
}

/*----------------------------------------------------------------------------*/

void IChainMatchExp::print(DebugOut &dbg) const {
   uint32_t ind = 2;
   dbg("chains assigned per class: ");
   dbg("%*s B1: %s", ind, "", chainsToStr(B1(), false).c_str());
   dbg("%*s B2: %s", ind, "", chainsToStr(B2()).c_str());
   dbg("%*s BX: %s", ind, "", chainsToStr(BX()).c_str());
   dbg("%*s BS: %s", ind, "", chainsToStr(BS()).c_str());
   dbg("%*s U1: %s", ind, "", chainsToStr(U1(), false).c_str());
   dbg("%*s UX: %s", ind, "", chainsToStr(UX()).c_str());
   dbg("%*s UC: %s", ind, "", chainsToStr(UC()).c_str());
   dbg("%*s US: %s", ind, "", chainsToStr(US()).c_str());
}

/*----------------------------------------------------------------------------*/

const Chains & IChainMatchExp::B1() const {
   return m_B1;
}

/*----------------------------------------------------------------------------*/

const Chains & IChainMatchExp::B2() const {
   return m_B2;
}

/*----------------------------------------------------------------------------*/

const Chains & IChainMatchExp::U1() const {
   return m_U1;
}

/*----------------------------------------------------------------------------*/

const Chains & IChainMatchExp::UX() const {
   return m_UX;
}

/*----------------------------------------------------------------------------*/

const Chains & IChainMatchExp::BX() const {
   return m_BX;
}

/*----------------------------------------------------------------------------*/

const Chains & IChainMatchExp::BS() const {
   return m_BS;
}

/*----------------------------------------------------------------------------*/

const Chains & IChainMatchExp::UC() const {
   return m_UC;
}

/*----------------------------------------------------------------------------*/

const Chains & IChainMatchExp::US() const {
   return m_US;
}
