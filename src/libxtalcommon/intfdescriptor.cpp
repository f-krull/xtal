#include "intfdescriptor.h"
#include "../libxtalcommon/chaininterfacetable.h"
#include "../xtalcompunbound/configcompunbound.h"
#include "../libxtaldata/protein.h"
#include "../libxtalutil/log.h"
#include "../libxtalutil/common.h"
#include "maxassignment.h"

#include <sstream>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------*/

UqChain::UqChain(char n, const std::string &s, const std::string &intfname) {
   seq = s;
   name = n;
   intfSize = 0;
   intf.resize(seq.size());
   /* extract interface */
   for (uint32_t i = 0; i < seq.size(); i++) {
      intf[i] = isupper(seq[i]) ? true : false;
   }
   /* count number of interface residues */
   for (uint32_t i = 0; i < seq.size(); i++) {
      intfSize += intf[i] == true ? 1 : 0;
   }
   this->intfName = intfname;
}

/*----------------------------------------------------------------------------*/

std::string UqChain::getIntStr() const {
   std::string ret(intf.size(), '#');
   for (uint32_t i = 0; i < intf.size(); i++) {
      ret[i] = intf[i] == true ? seq[i] : '_';
   }
   return ret;
}

/*----------------------------------------------------------------------------*/

UqEntry::UqEntry() {
   intname = "-";
   name = "-";
   resolution = 99;
   year = 0;
   numNonHohLig = 99;
}

/*----------------------------------------------------------------------------*/

/* favor small resolutions, greater years and greater names */
bool UqEntry::operator<(const UqEntry &e) const {
   if (e.numNonHohLig != numNonHohLig) {
      return e.numNonHohLig > numNonHohLig;
   }
   if (e.resolution != resolution) {
      return e.resolution > resolution;
   }
   if (e.year != year) {
      return e.year < year;
   }
   if (e.name != name) {
      return e.name < name;
   }
   return e.intname < intname;
}

/*----------------------------------------------------------------------------*/

void UqEntry::print() {
   Log::dbg("name       %s", name.c_str());
   Log::dbg("title      %s", title.c_str());
   Log::dbg("date       %02u", year);
   Log::dbg("resolution %4.2f", resolution);
   Log::dbg("found1     %s", chains1.c_str());
   Log::dbg("found2     %s", chains2.c_str());
   assert(chains1.size() == seq1.size());
   for (uint32_t i = 0; i < seq1.size(); i++) {
      Log::dbg("seq1       %c %s", chains1[i], seq1[i].seq.c_str());
      Log::dbg("int1       %c %s", chains1[i], seq1[i].getIntStr().c_str());
   }
   assert(chains2.size() == seq2.size());
   for (uint32_t i = 0; i < seq2.size(); i++) {
      Log::dbg("seq2       %c %s", chains2[i], seq2[i].seq.c_str());
      Log::dbg("int2       %c %s", chains2[i], seq2[i].getIntStr().c_str());
   }
}

/*----------------------------------------------------------------------------*/

static std::string chains2Str(const Chains &c) {
   std::string s;

   for (uint32_t i = 0; i < c.size(); i++) {
      s.push_back(c[i]->name());
   }
   return s;
}

/*----------------------------------------------------------------------------*/
#include <algorithm>
bool UqEntry::read(const Protein &pBC, const Chains &cB1, const Chains &cB2, const IntfDef &id) {
   /* compute interface */
   ChainInterfaceTable cIntfTab(id);
   cIntfTab.add(cB1, cB2);

   std::string chains1 = chains2Str(cB1);
   std::string chains2 = chains2Str(cB2);
   std::sort(chains1.begin(), chains1.end());
   std::sort(chains2.begin(), chains2.end());

   this->year = atoi(pBC.pdbInfo().date + 7);
   /* we introduce a year bug here - but its PDBs fault */
   this->year = this->year + (this->year > 50 ? 1900 : 2000);
   this->intname = name + " " + chains1 + " " + chains2;
   this->resolution = atof(pBC.pdbInfo().resolution);
   this->title = pBC.pdbInfo().title;

   for (uint32_t i = 0; i < cB1.size(); i++) {
      seq1.push_back(
            UqChain(cB1[i]->name(), cIntfTab.getSeqInt(cB1[i], cB2), intname));
      /* avoid zero length chains (caused by hetatm) */
      if (seq1.back().seq.size() <= 2) {
         seq1.pop_back();
      } else {
         this->chains1.push_back(cB1[i]->name());
      }
   }
   for (uint32_t i = 0; i < cB2.size(); i++) {
      seq2.push_back(
            UqChain(cB2[i]->name(), cIntfTab.getSeqInt(cB2[i], cB1), intname));
      /* avoid zero length chains */
      if (seq2.back().seq.size() <= 2) {
         seq2.pop_back();
      } else {
         this->chains2.push_back(cB2[i]->name());
      }
   }
   return this->chains1.empty() == false && this->chains2.empty() == false;
}

/*----------------------------------------------------------------------------*/

bool UqEntry::read(const std::string &fnpBC,
      const std::string &cnB1, const std::string &cnB2, const IntfDef &id) {
   Protein pBC;
   if (pBC.readPdb(fnpBC.c_str()) == false) {
      Log::err("error reading PDB file %s", fnpBC.c_str());
      return false;
   }

   ConfigCompUnbound cfg;

   /* select interfacing chains */
   Chains cB1 = IChains::getChainsByName(pBC, cnB1);
   Chains cB2 = IChains::getChainsByName(pBC, cnB2);

   /* check if selected */
   if (cB1.empty() == true || cB2.empty() == true) {
      Log::err("error %d", __LINE__);
      return false;
   }

   this->name = common::removeExtension(common::removePath(fnpBC));

   return read(pBC, cB1, cB2, id);
}

/*----------------------------------------------------------------------------*/

static uint32_t getTotalIntSize(const std::vector<UqChain> &seq1) {
   uint32_t numIntRes = 0;
   for (uint32_t i = 0; i < seq1.size(); i++) {
      numIntRes += seq1[i].intfSize;
   }
   return numIntRes;
}
/*----------------------------------------------------------------------------*/

static uint32_t getTotalSeqSize(const std::vector<UqChain> &seq1) {
   uint32_t numRes = 0;
   for (uint32_t i = 0; i < seq1.size(); i++) {
      numRes += seq1[i].seq.size();
   }
   return numRes;
}

/*----------------------------------------------------------------------------*/


#include "intfseqcomparer.h"

#define SCORE(x) x.getIntScore()

/*----------------------------------------------------------------------------*/


/* not sure if we should take min or max */
//#define GROUP_SMALL_WITH_BIG

#ifdef GROUP_SMALL_WITH_BIG
#define INT_CMP std::min
#else
#define INT_CMP std::max
#endif

/* store alignments, once computed */
class SimCache {
public:
   SimCache(const std::vector<UqChain> &s1, const std::vector<UqChain> &s2, bool ignoreSizeDiff);
   float getIntScore(const UqChain*, const UqChain*);
   std::map<std::pair<const UqChain*, const UqChain*>, float> m_mIntData;
   float getSeqScore(const UqChain*, const UqChain*);
   std::map<std::pair<const UqChain*, const UqChain*>, float> m_mSeqData;
private:

   uint32_t m_seqSize1;
   uint32_t m_seqSize2;
   uint32_t m_intSize1;
   uint32_t m_intSize2;
   bool m_ignoreSizeDiff;
};

SimCache::SimCache(const std::vector<UqChain> &seq1, const std::vector<UqChain> &seq2, bool ignoreSizeDiff) {

   /* count size of all interfaces */
   m_intSize1 = getTotalIntSize(seq1);
   m_seqSize1 = getTotalSeqSize(seq1);

   /* count size of all sequences */
   m_intSize2 = getTotalIntSize(seq2);
   m_seqSize2 = getTotalSeqSize(seq2);
   m_ignoreSizeDiff  = ignoreSizeDiff;
}

float SimCache::getIntScore(const UqChain* u1, const UqChain* u2) {
   std::map<std::pair<const UqChain*, const UqChain*>, float>::const_iterator it;
   it = m_mIntData.find(std::make_pair(u1, u2));
   if (it == m_mIntData.end()) {
      IntfSeqComparer isc;
      if (m_ignoreSizeDiff == false) {
         isc.align(u1->seq, u2->seq, std::max(m_seqSize1, m_seqSize2), std::max(m_intSize1, m_intSize2));
      } else {
         isc.align(u1->seq, u2->seq, std::min(m_seqSize1, m_seqSize2), std::min(m_intSize1, m_intSize2));
      }
      m_mIntData[std::make_pair(u1, u2)] = isc.getIntScore();
      m_mSeqData[std::make_pair(u1, u2)] = isc.getSeqId();
      return isc.getIntScore();
   }
   return it->second;
}

float SimCache::getSeqScore(const UqChain* u1, const UqChain* u2) {
   std::map<std::pair<const UqChain*, const UqChain*>, float>::const_iterator it;
   it = m_mSeqData.find(std::make_pair(u1, u2));
   if (it == m_mSeqData.end()) {
      IntfSeqComparer isc;
      if (m_ignoreSizeDiff == false) {
         isc.align(u1->seq, u2->seq, std::max(m_seqSize1, m_seqSize2), std::max(m_intSize1, m_intSize2));
      } else {
         isc.align(u1->seq, u2->seq, std::min(m_seqSize1, m_seqSize2), std::min(m_intSize1, m_intSize2));
      }
      m_mIntData[std::make_pair(u1, u2)] = isc.getIntScore();
      m_mSeqData[std::make_pair(u1, u2)] = isc.getSeqId();
      return isc.getIntScore();
   }
   return it->second;
}

/* comparators */
class CompUInt : public mass::CompU {
public:
   CompUInt(SimCache &m) : m_mc(m) {}

   float operator()(const UqChain* u1, const UqChain* u2) {
      return m_mc.getIntScore(u1, u2);
   }
private:
   SimCache &m_mc;
};
class CompUSeq : public mass::CompU {
public:
   CompUSeq(SimCache &m) : m_mc(m) {}

   float operator()(const UqChain* u1, const UqChain* u2) {
      return m_mc.getSeqScore(u1, u2);
   }
private:
   SimCache &m_mc;
};

/*----------------------------------------------------------------------------*/

class IntSeqMatch {
public:
   float getIntScore() const {return m_matchInt.score;}
   float getSeqId() const {return m_matchSeq.score;}

   void debug(const void *p) const;


   IntSeqMatch(const std::vector<UqChain> &seq1,
         const std::vector<UqChain> &seq2, bool ignoreSizeDiff);

private:
   mass::Match<UqChain> m_matchInt;
   mass::Match<UqChain> m_matchSeq;
};


void IntSeqMatch::debug(const void *p) const {
   for (uint32_t i = 0; i < m_matchInt.matches.size(); i++) {
      IntfSeqComparer isc;
      isc.label(m_matchInt.matches[i].s1->intfName, m_matchInt.matches[i].s1->intfName, m_matchInt.matches[i].s1->name, m_matchInt.matches[i].s1->name);
      isc.alignMin(m_matchInt.matches[i].s1->seq,m_matchInt.matches[i].s2->seq);
      isc.debug(p);
   }
}

/*----------------------------------------------------------------------------*/

IntSeqMatch::IntSeqMatch(const std::vector<UqChain> &seq1,
      const std::vector<UqChain> &seq2, bool ignoreSizeDiff) {

   SimCache mc(seq1, seq2, ignoreSizeDiff);

   const uint32_t lim = 20;
   CompUInt cmpI(mc);
   CompUSeq cmpS(mc);
   if (seq1.size() + seq2.size() <= lim) {
      mass::maximumWeightedMatching(seq1, seq2, cmpI, m_matchInt);
      mass::maximumWeightedMatching(seq1, seq2, cmpS, m_matchSeq);
   } else {
      mass::stableMarriageMatching(seq1, seq2, cmpI, m_matchInt);
      mass::stableMarriageMatching(seq1, seq2, cmpS, m_matchSeq);
   }
}

/*----------------------------------------------------------------------------*/

static bool isSimCandidate(const IntSeqMatch &int1,
      const IntSeqMatch &int2, const float threshold) {
   //return SCORE(int1) > threshold || SCORE(int2) > threshold;
   bool ret = false;
   ret = ret || int1.getSeqId() > threshold;
   ret = ret || int2.getSeqId() > threshold;
   ret = ret || int1.getIntScore() > threshold;
   ret = ret || int2.getIntScore() > threshold;
   return ret;
}

/*----------------------------------------------------------------------------*/

float DistUqentryCmp::operator ()(const UqEntry &e1, const UqEntry &e2) const {
   /* both sequences need to overcome threshold to be considered identical */
   const float SIM_THRESHOLD = 0.2;


   /* which match is best? calculate A->A + B->B */
   IntSeqMatch max1A(e1.seq1, e2.seq1, m_ignoreSizeDiff);
   IntSeqMatch max1B(e1.seq2, e2.seq2, m_ignoreSizeDiff);
   float simAB = 0;
   float simABseq = 0;
   float simABint = 0;
   simAB = IntfSeqComparer::combineIntSim(SCORE(max1A), SCORE(max1B));
   simABseq = IntfSeqComparer::combineIntSim(max1A.getSeqId(), max1B.getSeqId());
   simABint = IntfSeqComparer::combineIntSim(max1A.getIntScore(), max1B.getIntScore());

   /* is A->B + B->A better? */
   IntSeqMatch max2A(e1.seq1, e2.seq2, m_ignoreSizeDiff);
   IntSeqMatch max2B(e1.seq2, e2.seq1, m_ignoreSizeDiff);
   float simBA = 0;
   float simBAseq = 0;
   float simBAint = 0;
   simBA = IntfSeqComparer::combineIntSim(SCORE(max2A), SCORE(max2B));
   simBAseq = IntfSeqComparer::combineIntSim(max2A.getSeqId(), max2B.getSeqId());
   simBAint = IntfSeqComparer::combineIntSim(max2A.getIntScore(), max2B.getIntScore());

   /* pick better score */
   float sim = 0;
   enum {
      INT_AB,
      INT_BA
   } intOrder = INT_AB;

   if (simAB >= simBA) {
      if (m_showAlignment == true) {
         max1A.debug((void*)this);
         max1B.debug((void*)this);
      }
      intOrder = INT_AB;
      sim = simAB;
   } else {
      if (m_showAlignment == true) {
         max2A.debug((void*)this);
         max2B.debug((void*)this);
      }
      intOrder = INT_BA;
      sim = simBA;
   }

   /* only show if at least one score is not 0 */
   if (intOrder == INT_AB ?
         isSimCandidate(max1A, max1B, SIM_THRESHOLD) :
         isSimCandidate(max2A, max2B, SIM_THRESHOLD)) {
         Log::inf("intsc: %6.3f  seqsc: %6.3f  %s %s",
            intOrder == INT_AB ? simABint : simBAint,
            intOrder == INT_AB ? simABseq : simBAseq,
            (e1.name + "_" + e1.chains1 + "_" + e1.chains2).c_str(),
            (e2.name + "_" + e2.chains1 + "_" + e2.chains2).c_str()
            );
   }

   /* convert similarity to distance */
   return 1.0f - sim;
}

