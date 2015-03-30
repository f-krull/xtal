#include "cofactorgrouping.h"

#include <string.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <typeinfo>

#include "../libxtalutil/log.h"


#define MATCH_CAND_BY_LOWEST_GEOMETRIC_DIST 1

/*----------------------------------------------------------------------------*/

bool CofactorType::readGroupDef(const char* fn) {
   m_groups.clear();
   std::ifstream infile(fn);
   while (infile) {
      std::string s;
      if (getline(infile, s) == NULL) {
         break;
      }

      if (s.size() == 0 || s.at(0) == '#') {
         continue;
      }

      std::istringstream ss(s);
      std::set<std::string> record;

      while (ss) {
         std::string s;
         if (getline(ss, s, ' ') == NULL) {
            break;
         }
         record.insert(s);
      }
      if (record.empty() == false) {
         m_groups.push_back(record);
      }
   }
   if (infile.eof() == false) {
      infile.close();
      return false;
   }
   infile.close();

#if 0
   for (uint32_t i = 0; i < m_groups.size(); i++) {
      std::set<std::string>::iterator it = m_groups[i].begin();
      while (it != m_groups[i].end()) {
         printf("cof group %u : %s\n", i, it->c_str());
         it++;
      }
   }
#endif

   return m_groups.size() > 0;
}

CofactorType::MatchType CofactorType::matches(const Residue* c1, const Residue* c2) const {
   if (matchByName(c1, c2)) {
      return MATCH_BY_NAME;
   }
   if (matchByGroup(c1, c2)) {
      return MATCH_BY_GROUP;
   }
   if (matchByAtoms(c1, c2)) {
      return MATCH_BY_ATOMS;
   }
   return MATCH_NONE;
}

bool CofactorType::matchByName(const Residue* c1, const Residue* c2) const {
   return strcmp(c1->getName(), c2->getName()) == 0;
}

#define TO_STR(arg) #arg

const char* CofactorType::getMatchTypeStr(MatchType t) {
   switch (t) {
   case MATCH_BY_NAME:
      return TO_STR(MATCH_BY_NAME);
      break;
   case MATCH_BY_GROUP:
      return TO_STR(MATCH_BY_GROUP);
      break;
   case MATCH_BY_ATOMS:
      return TO_STR(MATCH_BY_ATOMS);
      break;
   default:
      return TO_STR(MATCH_NONE);
      break;
   }
   return TO_STR(MATCH_NONE);
}

#include "../libxtalutil/correlation.h"
#include <map>

static const char* trimFront(const char* s) {
   while (s[0] == ' ') {
      s++;
   }
   return s;
}

float getAtomCountSim(const Residue* c1, const Residue* c2) {

   std::map<char, std::pair<uint32_t, uint32_t> > h1h2;

   for (uint32_t i = 0; i < c1->atoms().size(); i++) {
      if (c1->atoms()[i]->name()[1] != 'H') {
         h1h2[c1->atoms()[i]->name()[1]].first++;
      }
   }
   for (uint32_t i = 0; i < c2->atoms().size(); i++) {
      if (c2->atoms()[i]->name()[1] != 'H') {
         h1h2[c2->atoms()[i]->name()[1]].second++;
      }
   }

   std::vector<uint32_t> hist1;
   std::vector<uint32_t> hist2;
   std::map<char, std::pair<uint32_t, uint32_t> >::const_iterator it;
   for (it = h1h2.begin(); it != h1h2.end(); ++it) {
      hist1.push_back(it->second.first);
      hist2.push_back(it->second.second);
   }

   int32_t count = 0; /* total number of atoms in a+b */
   int32_t diff = 0; /* total number of different atom counts */
   for (uint32_t i = 0; i < hist1.size(); i++) {
      count += hist1[i] + hist2[i];
      diff += abs(hist1[i] - hist2[i]);
   }
   float similarity = 1 - (float(diff) / count);

   /* penalize high overall differences in number of atoms */
   similarity *= std::min(c1->atoms().size(), c2->atoms().size());
   similarity /= std::max(c1->atoms().size(), c2->atoms().size());

#if 0

   printf("correlation %s_%s %f\n", trimFront(c1->getName()), trimFront(c2->getName()),
         similarity);
   for (it = h1h2.begin(); it != h1h2.end(); ++it) {
      printf("%c:%u/%u ", it->first, it->second.first, it->second.second);
   }

   printf("correlation %s %s %f\n", c1->getName(), c2->getName(),
         corrSpears(hist1, hist2));

   for (uint32_t i = 0; i < hist1.size(); i++) {
      printf("%2.0f ", hist1[i]);
   }
   printf("\n");
   for (uint32_t i = 0; i < hist2.size(); i++) {
      printf("%2.0f ", hist2[i]);
   }
   printf("\n");
#endif
   return similarity;
}

bool CofactorType::matchByAtoms(const Residue* c1, const Residue* c2) const {

   /* correlation of the two histograms */
   return (getAtomCountSim(c1, c2) > 0.6);
}

static void remSpaces(std::string &s) {
   s.erase(std::remove_if(s.begin(), s.end(), isspace), s.end());
}

bool CofactorType::matchByGroup(const Residue* c1, const Residue* c2) const {
   std::string k1 = c1->getName();
   remSpaces(k1);
   for (uint32_t i = 0; i < m_groups.size(); i++) {
      /* search first cofactor name in current group */
      if (m_groups[i].find(k1) == m_groups[i].end()) {
         continue;
      }
      /* search cofactor name */
      std::string k2 = c2->getName();
      remSpaces(k2);
      if (m_groups[i].find(k2) == m_groups[i].end()) {
         continue;
      }
      /* match */
      return true;
   }

   return false;
}

/*----------------------------------------------------------------------------*/

#include <set>

#include "../libxtalutil/common.h"

CofFamily::CofFamily(Residue * r) {
   add(r);
}

bool CofFamily::matches(const Residue* r, const CofactorType &c) const {
   bool found = false;
   /* see if it matches initial cofactor */
   found = found || c.matches(m_cofactors.front(), r);
   return found;
}

void CofFamily::add(Residue* r) {
   m_types.insert(r->getName());
   m_cofactors.push_back(r);
}

std::string CofFamily::type() const {
   return trimFront(m_cofactors.front()->getName());
}

std::string CofFamily::toStr() const {
   std::string ret;
   ret += trimFront(m_cofactors.front()->getName());
   ret += "(";
   ret += common::s_printf("%u:", m_cofactors.size());
   for (std::set<CofType>::const_iterator it = m_types.begin();
         it != m_types.end(); ++it) {
      ret += (it == m_types.begin()) ? "" : ",";
      ret += trimFront(it->c_str());
   }
   ret += ")";
   return ret;
}

const Residues & CofFamily::cofactors() const {
   return m_cofactors;
}

/*----------------------------------------------------------------------------*/

static bool cmpSmaller(const Residue* a, const Residue* b) {
   return strcmp(a->getName(), b->getName()) == -1;
}

CofGrouping::CofGrouping(const CofactorType &cm, const Residues &cof) {
   Residues c(cof.begin(), cof.end());
   std::sort(c.begin(), c.end(), cmpSmaller);
   for (uint32_t i = 0; i < c.size(); i++) {
      bool found = false;
      /* try to add to existing family */
      for (uint32_t j = 0; j < m_families.size(); j++) {
         if (m_families[j].matches(c[i], cm) == true) {
            m_families[j].add(c[i]);
            found = true;
            break;
         }
      }
      /* if no family matches, create own */
      if (found == false) {
         m_families.push_back(c[i]);
      }
   }
}

std::string CofGrouping::toStr() const {
   std::string ret = ";";
   for (uint32_t i = 0; i < m_families.size(); i++) {
      ret += m_families[i].toStr();
      ret += ";";
   }
   return ret;
}

const std::vector<CofFamily> & CofGrouping::families() const {
   return m_families;
}

/*----------------------------------------------------------------------------*/

#include <assert.h>

typedef char cbool;

static uint32_t matchCof(const CofactorDetector::CofNeighborhood &cofB1,
      const CofactorDetector::CofNeighborhood &cofU1,
      const std::map<Residue*, Residue*> &alignB1U1) {
   uint32_t numMatched = 0;

   /* check all neighbors of cofB1 */
   std::set<Residue*>::const_iterator it;
   for (it = cofB1.neighbors.begin(); it != cofB1.neighbors.end(); ++it) {
      /* try find the aligned U1 residue */
      std::map<Residue*, Residue*>::const_iterator mU1 = alignB1U1.find(*it);
      if (mU1 == alignB1U1.end()) {
         continue;
      }
      /* try to find translated residue within U1 neighbors */
      if (cofU1.neighbors.find(mU1->second) == cofU1.neighbors.end()) {
         continue;
      }
      numMatched++;
   }
   return numMatched;
}

#include <queue>

struct MatchCand {
   uint32_t numMatched;
   uint32_t numTotal;
   uint32_t idxB;
   uint32_t idxU;
   Residue *resB;
   Residue *resU;
   float geomDist;
   float nuDivNb;
};



class MatchCandRevComp {
public:
   bool operator()(const MatchCand &a, const MatchCand &b) const {
#if MATCH_CAND_BY_LOWEST_GEOMETRIC_DIST
      if (a.geomDist / a.nuDivNb  != b.geomDist / b.nuDivNb) {
         return a.geomDist / a.nuDivNb > b.geomDist / b.nuDivNb;
      }
#endif
      if (a.numMatched != b.numMatched) {
         return a.numMatched < b.numMatched;
      }
      if (a.resB->atoms().front()->chainId()
            != b.resB->atoms().front()->chainId()) {
         return a.resB->atoms().front()->chainId()
               < b.resB->atoms().front()->chainId();
      }
      if (strcmp(a.resB->atoms().front()->resName(),
            b.resB->atoms().front()->resName()) != 0) {
         return strcmp(a.resB->atoms().front()->resName(),
               b.resB->atoms().front()->resName()) < 0;
      }
      if (a.resB->atoms().front()->iCode()
            != b.resB->atoms().front()->iCode()) {
         return a.resB->atoms().front()->iCode()
               < b.resB->atoms().front()->iCode();
      }
      if (a.resU->atoms().front()->chainId()
            != b.resU->atoms().front()->chainId()) {
         return a.resU->atoms().front()->chainId()
               < b.resU->atoms().front()->chainId();
      }
      if (strcmp(a.resU->atoms().front()->resName(),
            b.resU->atoms().front()->resName()) != 0) {
         return strcmp(a.resU->atoms().front()->resName(),
               b.resU->atoms().front()->resName()) < 0;
      }
      if (a.resU->atoms().front()->iCode()
            != b.resU->atoms().front()->iCode()) {
         return a.resU->atoms().front()->iCode()
               < b.resU->atoms().front()->iCode();
      }

      return a.resB < b.resB;
   }
};


static Vector resCentre(const Residue* r) {
   Vector centre(0, 0, 0);
   unsigned int num;

   num = 0;
   for (uint32_t j = 0; j < r->atoms().size(); j++) {
      centre += (*r->atoms()[j]);
      num++;
   }
   if (num != 0) {
      centre /= num;
      return centre;
   }
   return Vector(0, 0, 0);
}

static float getGeomDist(const Residue *c1, const Residue *c2) {
   return resCentre(c1).dist(resCentre(c2));
}

bool isValidMatch(const MatchCand &c, uint32_t minNumMatched) {
   /* anything matched? */
   if (c.numMatched == 0) {
      return false;
   }
   /* enough matched? */
   if (c.numMatched < minNumMatched) {
      return false;
   }
   return true;
}


CofactorMatcher::CofactorMatcher(const CofactorType &c,
      const std::vector<CofactorDetector::CofNeighborhood> &cofB,
      const Residues &intB,
      const std::vector<CofactorDetector::CofNeighborhood> &cofU,
      const Residues &intU, const CofactorType &ct, DebugOut &dbg) {
   m_numMatched = 0;
   m_numUnmatched = cofB.size();

   /* create map to find B1 residues faster */
   std::map<Residue*, Residue*> alignB1U1;
   assert(intB.size() == intU.size());
   for (uint32_t i = 0; i < intB.size(); i++) {
      alignB1U1[intB[i]] = intU[i];
   }

   std::priority_queue<MatchCand, std::vector<MatchCand>, MatchCandRevComp> matchCands;

   dbg("cofactor cand sizes %u(b) -> %u(u)", cofB.size(),
         cofU.size());
   for (uint32_t i = 0; i < cofB.size(); i++) {
      for (uint32_t j = 0; j < cofU.size(); j++) {
         MatchCand mc;
         mc.idxB = i;
         mc.idxU = j;
         mc.resB = cofB[i].cof;
         mc.resU = cofU[j].cof;
         mc.numTotal = cofB[i].neighbors.size();
         mc.geomDist = getGeomDist(mc.resB, mc.resU);
         /* has similar neighbors? */
         mc.numMatched = matchCof(cofB[i], cofU[j], alignB1U1);
         mc.nuDivNb = float(mc.numMatched) / mc.numTotal;
         uint32_t minNumMatched = std::max((uint32_t)cofB[i].neighbors.size() / 4u, 1u);
         if (mc.numMatched < minNumMatched) {
            dbg("match candidate %s -> %s: num.Neigh.B=%u num.Neigh.U=%u (%4.2f)", mc.resB->getName(), mc.resU->getName(), mc.numTotal, mc.numMatched, float (mc.numMatched)/mc.numTotal);
            continue;
         }
         const CofactorType::MatchType matches = ct.matches(mc.resB, mc.resU);
         if (matches != CofactorType::MATCH_NONE) {
            matchCands.push(mc);
         }
         dbg("match candidate %s -> %s: %s Nb=%2u Nb/Nu=%4.2f dist=%.1f", mc.resB->getName(), mc.resU->getName(), CofactorType::getMatchTypeStr(matches), mc.numTotal, float (mc.numMatched)/mc.numTotal, mc.geomDist);
      }
   }

   /* do not reassign U1 residues */
   std::vector<cbool> cofBmatched(cofB.size(), false);
   std::vector<cbool> cofUmatched(cofU.size(), false);

   assert(cofUmatched.size() == cofU.size());
   assert(cofUmatched.size() == cofU.size());
   /* try to add candidates */
   while (matchCands.empty() == false) {
      /* already matched? */
      if (cofBmatched[matchCands.top().idxB] == false
            && cofUmatched[matchCands.top().idxU] == false) {
         /* add */
         m_cofMatch[matchCands.top().resB] = matchCands.top().resU;
         m_cofMatch[matchCands.top().resU] = matchCands.top().resB;
         cofUmatched[matchCands.top().idxU] = true;
         cofBmatched[matchCands.top().idxB] = true;
         dbg("found match %s to %s with %u of %u sim: %5.2f dist:%5.1f", matchCands.top().resB->getName(), matchCands.top().resU->getName(), matchCands.top().numMatched, matchCands.top().numTotal, getAtomCountSim(matchCands.top().resB, matchCands.top().resU), matchCands.top().geomDist);
         m_numMatched++;
         m_numUnmatched--;
      }
      matchCands.pop();
   }

   assert(m_numMatched + m_numUnmatched == cofB.size());
}

uint32_t CofactorMatcher::numMatched() const {
   return m_numMatched;
}
uint32_t CofactorMatcher::numUnmatched() const {
   return m_numUnmatched;
}

Residue* CofactorMatcher::getMatched(Residue *res) const {
   std::map<Residue*, Residue*>::const_iterator it = m_cofMatch.find(res);
   if (it != m_cofMatch.end()) {
      return it->second;
   }
   return NULL;
}
