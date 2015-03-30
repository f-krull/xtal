#include "cofactoruniqdescriptor.h"

/*----------------------------------------------------------------------------*/

static const char* trimFront(const char* s) {
   while (s[0] == ' ') {
      s++;
   }
   return s;
}

/*----------------------------------------------------------------------------*/

static std::string resToSTr(const Residue* r) {
   std::string s = "";
   s += r->atoms().front()->chainId();
   s += ":";
   s += trimFront(r->atoms().front()->resSeq());
   if (r->atoms().front()->iCode() != ' ') {
      s += r->atoms().front()->iCode();
   }
   s += std::string("(") + trimFront(r->atoms().front()->resName()) + ")";
   return s;
}

static void addCof(std::map<std::string, std::string> &sCofStr
      , const CofactorDetector &cd, const CofactorMatcher *match) {
   for (uint32_t i = 0; i < cd.cofactors().size(); i++) {
      std::string matchStr;
      if (match != NULL && match->getMatched(cd.cofactors()[i].cof) != NULL) {
         matchStr = resToSTr(match->getMatched(cd.cofactors()[i].cof));
      } else {
         matchStr = "";
      }
      sCofStr[resToSTr(cd.cofactors()[i].cof)] = matchStr;
   }
}

UniqueCofactorDescriptor::UniqueCofactorDescriptor(const CofactorDetector &dB1, const CofactorDetector &dB2, const CofactorDetector &dU1, const CofactorMatcher &match) {

   /* unique residue - custom status */
   std::map<std::string, std::string> sCofStr;

#if 0
   if (m_cofactors.empty() == false) {
      ret += m_cofactors.front()->atoms().front()->chainId();
      ret += ":";
      ret += trimFront(m_cofactors.front()->atoms().front()->resSeq());
      if (m_cofactors.front()->atoms().front()->iCode() != ' ') {
         ret += m_cofactors.front()->atoms().front()->iCode();
      }
   }
#endif

   addCof(sCofStr, dB2, NULL); /* B2 cofactors are not matched! */
   addCof(sCofStr, dB1, &match); /* allow overwrite of entries with status */

   m_cofStr= ";";
   std::map<std::string, std::string>::const_iterator it;
   for (it = sCofStr.begin(); it != sCofStr.end(); it++) {
      m_cofStr += it->first + "," + it->second + ";";
   }

   /* add unmatched U1 cofactors as well */
   for (uint32_t i = 0; i < dU1.cofactors().size(); ++i) {
      if (match.getMatched(dU1.cofactors()[i].cof) == NULL) {
         m_cofStr += "," + resToSTr(dU1.cofactors()[i].cof) + ";";
      }
   }

   m_numCofactors = sCofStr.size();
}

const std::string & UniqueCofactorDescriptor::cofStr() const {return m_cofStr;}

uint32_t UniqueCofactorDescriptor::numCofactors() const {return m_numCofactors;}
