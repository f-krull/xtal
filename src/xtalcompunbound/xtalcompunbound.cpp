#include "xtalcompunbound.h"
#include "debugoutput.h"
#include "../libxtalutil/command.h"
#include "../libxtalutil/common.h"
#include "../libxtalutil/chunkstatus.h"
#include "../libxtaldata/protein.h"
#include "chainmatch.h"
#include "multichaininterface.h"
#include "chainmatchexpansion.h"
#include "chainoverlapinfo.h"
#include "cofactordetector.h"
#include "cofactorgrouping.h"
#include "cofactoruniqdescriptor.h"

#include "resaligner.h"

#include "../libxtalcommon/interface.h"
#include "../libxtalcommon/seqaligner.h"
#include "../libxtalcommon/kabschwrapper.h"
#include "../libxtalcommon/chaininterfacetable.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <algorithm>
#include <memory>

#define ENV_COF_IGNORELIST "XTAL_COF_IGNORELIST"
#define ENV_COF_GROUPS "XTAL_COF_GROUPS"
#define MAXSEQID 0.7

/*----------------------------------------------------------------------------*/

const char* XtalCompUnbound::getName() const {
	return "XtalCompUnbound";
}

/*----------------------------------------------------------------------------*/

class ChainNames {
public:

   const std::vector<std::string> &b1() const;
   const std::vector<std::string> &u1() const;
   const std::vector<std::string> &b2() const;

   bool parse(std::string chainsb1, std::string chainsu1, std::string chainsb2);
private:
   std::vector<std::string> m_bound1;
   std::vector<std::string> m_unboud1;
   std::vector<std::string> m_bound2;
};

/*----------------------------------------------------------------------------*/

class UnboundMatchResult {
public:
   enum Status {
      ERR_OK = 0,
      ERR_CM_BAD = 1,
      ERR_CM_UC = 2,
      ERR_IMSD = 4,
      ERR_INTALGN = 8,
      ERR_DNA = 16,
      ERR_U1B2CLASH = 32,
      ERR_U1INTSIM = 64,
      ERR_NO_B1B2INT = 128,
      ERR_PDBFILE = 256,
      ERR_CHAINSB1 = 512,
      ERR_CHAINSB2 = 1024,
      ERR_CHAINSU1 = 2048,
      ERR_CHAINCLASHES = 4096,
      ERR_DNAU1 = 8192,
      ERR_BACKBONEB = 16384,
      ERR_BACKBONEU1 = 32768
   };

   uint32_t status;

   std::string seedPBC;
   std::string seedB1;
   std::string seedB2;
   std::string seedPU1;
   std::string seedU1;

   uint32_t cNumInterfaceChains;
   uint32_t cNumGaps;
   uint32_t cNumNonInterfaceChains;
   uint32_t cNumLigands;
   uint32_t cNumS2Bonds;

   std::string cChainsB1;
   std::string cChainsB2x;
   std::string cChainsB1i;
   std::string cChainsB2i;

   uint32_t cNumI1Gaps;
   uint32_t cNumI2Gaps;
   uint32_t cNumI1Ca;
   uint32_t cNumI2Ca;



   uint32_t u1NumChains;
   uint32_t u1NumAddChains;
   uint32_t u1NumGaps;
   float uAlignedIntRatio;
   uint32_t uNumMatchedLigands;
   uint32_t uNumUnmatchedLigands;
   float uIRmsd;
   uint32_t uNumClashes;

   std::string uChainsU1s;
   std::string uChainsU1x;

   std::vector<float> rotTrans;

   std::string cofStr;

   std::string errInfo;

   UnboundMatchResult() {
      cNumInterfaceChains = 0;
      cNumGaps = 0;
      cNumNonInterfaceChains = 0;
      cNumLigands = 0;
      cNumS2Bonds = 0;

      u1NumChains = 0;
      u1NumAddChains = 0;
      u1NumGaps = 0;
      uAlignedIntRatio = 0;
      uNumMatchedLigands = 0;
      uNumUnmatchedLigands = 0;
      uIRmsd = 0;
      uNumClashes = 0;

      cofStr = ";";

      cNumI1Gaps = 0;
      cNumI2Gaps = 0;
      cNumI1Ca = 0;
      cNumI2Ca = 0;

      rotTrans.resize(12, 0.0f);
      status = ERR_OK;
   }


   std::string toTableRow() const {
      std::string res;
      res += common::s_printf("%5s ", seedPBC.c_str());
      res += common::s_printf("%s ", seedB1.c_str());
      res += common::s_printf("%s ", seedB2.c_str());
      res += common::s_printf("%5s ", seedPU1.c_str());
      res += common::s_printf("%s ", seedU1.c_str());

      res += common::s_printf("%s ", status == ERR_OK ? "ok" : "error");



      res += common::s_printf("%4s ", cChainsB1i.c_str());
      res += common::s_printf("%4s ", cChainsB2i.c_str());

      res += common::s_printf("%2u ", cNumInterfaceChains);
      res += common::s_printf("%2u ", cNumGaps);
      res += common::s_printf("%2u ", cNumI1Gaps);
      res += common::s_printf("%2u ", cNumI2Gaps);

      res += common::s_printf("%2u ", cNumI1Ca);
      res += common::s_printf("%2u ", cNumI2Ca);


      res += common::s_printf("%2u ", cNumNonInterfaceChains);
      res += common::s_printf("%2u ", cNumLigands);
      res += common::s_printf("%2u ", cNumS2Bonds);

      res += common::s_printf("%6s ", cChainsB1.c_str());
      res += common::s_printf("%6s ", cChainsB2x.c_str());

      res += common::s_printf("%2u ", u1NumChains);
      res += common::s_printf("%2u ", u1NumGaps);
      res += common::s_printf("%2u ", u1NumAddChains);

      res += common::s_printf("%4.2f ", uAlignedIntRatio);

      res += common::s_printf("%2u ", uNumMatchedLigands);
      res += common::s_printf("%2u ", uNumUnmatchedLigands);

      res += common::s_printf("%4.2f ", uIRmsd);
      res += common::s_printf("%2u ", uNumClashes);

      res += common::s_printf("%4s ", uChainsU1x.c_str());

      for (uint32_t i = 0; i < rotTrans.size(); i++) {
         res += common::s_printf("%6.2f ", rotTrans[i]);
      }

      res += common::s_printf("%s ", cofStr.c_str());

      std::string errStr;
      if (status & ERR_CHAINCLASHES) {
         errStr += common::s_printf("ERR_CHAINCLASHES - one PDB structure has overlapping chains");
      } else if (status & ERR_NO_B1B2INT) {
         errStr += common::s_printf("ERR_NO_B1B2INT B1:B2 (seed) do not have an interface.");
      } else if (status & ERR_CM_BAD) {
         errStr += common::s_printf("ERR_CM_BAD - no chainmatch was found at all.");
      } else if (status & ERR_CM_UC) {
         errStr += common::s_printf("ERR_CM_UC - some U1 chains are not covered by chainmatch, but have an interface to B2.");
      } else if (status & ERR_DNA) {
         errStr += common::s_printf("ERR_DNA - DNA was found within interface of PB.");
      } else if (status & ERR_INTALGN) {
         errStr += common::s_printf("ERR_INTALGN - U1 interface residues do not match well to B1.");
      } else if (status & ERR_IMSD) {
         errStr += common::s_printf("ERR_IMSD - high iRMSD");
      } else if (status & ERR_U1B2CLASH) {
         errStr += common::s_printf("ERR_U1B2CLASH - U1 chains clash with B2 after superimposition.");
      } else if (status & ERR_U1INTSIM) {
         errStr += common::s_printf("ERR_U1INTSIM - an interface similar to B1:B2 (seed) was found in PU.");
      } else if (status & ERR_PDBFILE) {
         errStr += common::s_printf("ERR_PDBFILE - PDB file was not found. Something went very wrong! ");
      } else if (status & ERR_DNAU1) {
         errStr += common::s_printf("ERR_DNAU1 - DNA was found within interface of PU.");
      } else if (status & ERR_BACKBONEB) {
         errStr += common::s_printf("ERR_BACKBONEB - bound PDB structure is missing many backbone atoms");
      } else if (status & ERR_BACKBONEU1) {
         errStr += common::s_printf("ERR_BACKBONEU1 - unbound PDB structure is missing many backbone atoms");
      } else if (status != ERR_OK) {
         errStr += common::s_printf("ERR_%d(not handled) ", status);
      }

      if (status != ERR_OK) {
         if (errInfo.empty() == false) {
            errStr += common::s_printf("(%s)", errInfo.c_str());
         }
      }
      std::replace(errStr.begin(), errStr.end(), ' ', '_');
      res += errStr;


      return res;
   }
};

/*----------------------------------------------------------------------------*/

static std::vector<std::string> tokenize(const std::string &inpStr, char delim =
         ' ') {
   std::vector<std::string> args;
   std::string token;
   std::istringstream iss(inpStr);
   while (std::getline(iss, token, delim)) {
      if (token.size() != 0) {
         args.push_back(token);
      }
   }
   return args;
}

/*----------------------------------------------------------------------------*/

const std::vector<std::string> & ChainNames::b1() const {
   return m_bound1;
}

/*----------------------------------------------------------------------------*/

const std::vector<std::string> & ChainNames::u1() const {
   return m_unboud1;
}

/*----------------------------------------------------------------------------*/

const std::vector<std::string> & ChainNames::b2() const {
   return m_bound2;
}

/*----------------------------------------------------------------------------*/

static void removeChar(std::string &str, char c) {
   std::remove(str.begin(), str.end(), c);
}

/*----------------------------------------------------------------------------*/

bool ChainNames::parse(std::string chainsb1, std::string chainsb2,
         std::string chainsu1) {
   removeChar(chainsb1, '"');
   removeChar(chainsu1,'"');
   removeChar(chainsb2,'"');
   m_bound1 = tokenize(chainsb1, ',');
   m_unboud1 = tokenize(chainsu1, ',');
   m_bound2 = tokenize(chainsb2, ',');

   if (m_bound1.size() == 0 || m_unboud1.size() == 0 || m_bound2.size() == 0) {
      return false;
   }

   /* every chain got a translation */
   if (m_bound1.size() != m_unboud1.size()) {
      return false;
   }
   return true;
}

/*----------------------------------------------------------------------------*/

static bool isAlignment(const std::map<Residue*, bool> &m, Residue* r) {
   std::map<Residue*, bool>::const_iterator it = m.find(r);
   return it != m.end();
}

/*----------------------------------------------------------------------------*/
#include "../libxtaldata/aamap.h"
static void printIntAlignment(const ResAligner::Alignment &al, const Residues &resB1,
      const std::vector<cbool> intB1, DebugOut &dbg) {
   assert(resB1.size() == intB1.size());

   std::map<Residue*, bool> isIntMap;
   for (uint32_t i = 0; i < resB1.size(); ++i) {
      if (intB1[i] == true) {
         isIntMap[resB1[i]] = intB1[i];
      }
   }

   std::string sal1;
   std::string sal2;
   std::string sali;
   std::string saln;

   for (uint32_t i = 0; i < al.alfull1.size(); ++i) {
      const Residue *r = al.alfull1[i];
      sal1.push_back(r != NULL ? AaMap::getResn1(r->atoms()[0]->resName()) : '-');
   }
   for (uint32_t i = 0; i < al.alfull2.size(); ++i) {
      const Residue *r = al.alfull2[i];
      sal2.push_back(r != NULL ? AaMap::getResn1(r->atoms()[0]->resName()) : '-');
   }
   for (uint32_t i = 0; i < al.alfull1.size(); ++i) {
      sali.push_back(isAlignment(isIntMap, al.alfull1[i]) ? '*' : ' ');
   }
   for (uint32_t i = 0; i < al.match.size(); ++i) {
      saln.push_back(al.match[i] == true ? '*' : ' ');
   }
   dbg("intaln   \"int\": \"%s\",", sali.c_str());
   dbg("intaln   \"alb\": \"%s\",", sal1.c_str());
   dbg("intaln   \"alu\": \"%s\",", sal2.c_str());
   dbg("intaln   \"aln\": \"%s\"", saln.c_str());
}

/*----------------------------------------------------------------------------*/

static float getAlignedInterface(const Chains &chainsB1, const Chains &chainsU1,
      const std::vector<std::vector<cbool> > &intB1,
      std::vector<Residue*> &intResB1, std::vector<Residue*> &intResU1, DebugOut &dbg, const UnboundMatchResult &ur) {
	uint32_t numalint;
	uint32_t numint;
	uint32_t k;
	Residue* res;

   uint32_t numalint_total = 0;
   uint32_t numint_total = 0;

   bool isFirst = true;

   dbg("intaln {");

	/* align all chains */
	assert(chainsB1.size() == intB1.size());
	assert(chainsB1.size() == chainsU1.size());
	for (uint32_t i = 0; i < chainsB1.size(); ++i) {

      /* align all residues */
      ResAligner::Alignment al = ResAligner::align(chainsB1[i]->resis(), chainsU1[i]->resis());
#if 0
      /* debug */
      {
         std::string alStr1;
         std::string alStr2;

         for (uint32_t k = 0; k < al.al1.size(); k++) {
            alStr1.push_back(AaMap::getResn1(al.al1[k]->atoms()[0]->resName()));
            alStr2.push_back(AaMap::getResn1(al.al2[k]->atoms()[0]->resName()));
         }
         dbg("B1-%c: %s", chainsB1[i]->name(), alStr1.c_str());
         dbg("U1-%c: %s", chainsU1[i]->name(), alStr2.c_str());
      }
#endif

      /* count aligned interface residues */
      numalint = 0;
      numint = 0;
      k = 0;

      /* look at reach residue */
      for (uint32_t l = 0; l < chainsB1[i]->resis().size(); l++) {
         if (intB1[i][l] == true) {
            numint++;
         }
         if (k >= al.al1.size()) {
            continue;
         }
         res = chainsB1[i]->resis()[l];
         /* residue belongs to interface? */
         if (intB1[i][l] == true) {
            /* next aligned residue found? */
            if (res == al.al1[k]) {
               intResB1.push_back(al.al1[k]);
               intResU1.push_back(al.al2[k]);
               numalint++;
            }
         }
         /* found aligned residue -> go to next one */
         if (res == al.al1[k]) {
            k++;
         }
      }

      if (numint > 0) {
         dbg("intaln %c\"%d\" : {", isFirst == false ? ',' :  ' ', i);
         dbg("intaln   \"pdbb\": \"%s\",", ur.seedPBC.substr(0, 4).c_str());
         dbg("intaln   \"pdbu\": \"%s\",", ur.seedPU1.substr(0, 4).c_str());
         dbg("intaln   \"chainb\": \"%c\",", chainsB1[i]->name());
         dbg("intaln   \"chainu\": \"%c\",", chainsU1[i]->name());
         dbg("intaln   \"nint\": \"%d\",", numint);
         dbg("intaln   \"nalint\": \"%d\",", numalint);
         printIntAlignment(al, chainsB1[i]->resis(), intB1[i], dbg);
         dbg("intaln  }");
         isFirst = false;
      }


      dbg("chain %c -> %c : %u / %u bound interface residues covered by alignment with unbound", chainsB1[i]->name(), chainsU1[i]->name(),
            numalint, numint);

      numalint_total += numalint;
      numint_total += numint;
   }

	dbg("intaln }");

	return numint_total > 0 ? (float) numalint_total / numint_total : 0.0f;
}


/*----------------------------------------------------------------------------*/

static void getCaAtoms(const std::vector<Residue*> &rB1, std::vector<Atom*> &caB1, const std::vector<Residue*> &rU1, std::vector<Atom*> &caU1) {
   caB1.clear();
   caU1.clear();
   assert(rB1.size() == rU1.size());
   for (unsigned int i = 0; i < rB1.size(); i++) {
      if (rB1[i]->hasCa() == false || rU1[i]->hasCa() == false) {
         continue;
      }
      caB1.push_back(rB1[i]->ca());
      caU1.push_back(rU1[i]->ca());
   }
}

/*----------------------------------------------------------------------------*/
/* calculate the mean of all added vectors */
class VecCentre : public Vector {
public:
   VecCentre() {
      m_numData = 0;
   }

   template<typename T>
   void add(const std::vector<T> &data) {
      for (uint32_t i = 0; i < data.size(); i++) {
         m_sumData += (*data[i]);
         m_numData++;
      }
      Vector::operator=(m_sumData / m_numData);
   }

private:
   uint32_t m_numData;
   Vector m_sumData;
};

/*----------------------------------------------------------------------------*/

static std::string chains2Str(const Chains &c) {
   std::string s;

   for (uint32_t i = 0; i < c.size(); i++) {
      s.push_back(c[i]->name());
   }
   return s;
}

/*----------------------------------------------------------------------------*/

static std::string chains2Str_sorted(const Chains &c) {
   std::string s;

   for (uint32_t i = 0; i < c.size(); i++) {
      s.push_back(c[i]->name());
   }
   std::sort(s.begin(), s.end());
   return s;
}

/*----------------------------------------------------------------------------*/

#if 0
static bool hasChainSeqConfilct(const Protein &pBound,
         const std::vector<std::string> &cBound2, const Protein &pUnbound,
         const std::vector<std::string> &cUnbound1) {
   std::vector<Chain*> selBound;
   std::vector<Chain*> selUnbound;

   selBound= getChainsByName(pBound, cBound2);
   selUnbound = getChainsByName(pUnbound, cUnbound1);

   assert(selBound.size());
   assert(selUnbound.size());

   for (uint32_t i = 0; i < selBound.size(); i++) {
      std::string seqBound = Residue::getResSequence(selBound[i]->resis());
      /* happens if chain has no amino acid residues */
      if (seqBound.size() == 0) {
         continue;
      }
      for (uint32_t j = 0; j < selUnbound.size(); j++) {
         std::string seqUnbound = Residue::getResSequence(
                  selUnbound[i]->resis());
         /* happens if chain has no amino acid residues */
         if (seqUnbound.size() == 0) {
            continue;
         }

         /* is any Bound2 chain part of any Unbound chain? */
         float seqid = getSeqIdDirected(seqBound, seqUnbound);
         if (seqid >= MAXSEQID) {
         Log::dbg("seqid bound_forbidden %2u with unbound %2u: %5.2f chainlengths: %lu, %lu", i ,j, seqid, seqBound.size(),
                  seqUnbound.size());

            return true;
         }
      }
   }
   return false;
}
#endif

/*----------------------------------------------------------------------------*/

bool hasChainOverlap(const Chains &cB1, const Chain* cB2, const Chains &cU1, ChainOverlapInfo &ovlInfo, float u1ToB2MaxConflictRatio, DebugOut &dbg) {

   /* require correct alignment */
   assert(cU1.size() == cB1.size());

   for (uint32_t i = 0; i < cU1.size(); i++) {
      const float ovl = ovlInfo.getChainOverlap(cB1[i], cB2, cU1[i]);
      dbg("%s: %c->%c: %f", "hasChainOverlap", cU1[i]->name(), cB2->name(), ovl);
      if (ovl > u1ToB2MaxConflictRatio) {
         dbg("%s: clash! %c->%c: %f", "hasChainOverlap", cU1[i]->name(), cB2->name(), ovl);
         return true;
      }
   }

   return false;
}

/*----------------------------------------------------------------------------*/
#include "../libxtalcommon/intfdescriptor.h"
bool hasIntSimilarity(const Protein &pBC, const Chains &cB1, const Chains &cB2, const Protein &pU1, const Chains &cU1, float u1MinIntDist, bool debug = false) {

   Chains cU2; /* get all unbound chains that are not cU1 */
   {
      Chains cU1s = cU1;
      std::sort(cU1s.begin(), cU1s.end());
      for (uint32_t i = 0; i < pU1.chains().size(); i++) {
         /* only add if not found in U1 */
         if (std::binary_search(cU1s.begin(), cU1s.end(), pU1.chains()[i]) == false) {
            cU2.push_back(pU1.chains()[i]);
         }
      }
   }

   /* get interface similarity of B1:B2 with U1:U2 */
   UqEntry intB;
   UqEntry intU;
   if (intB.read(pBC, cB1, cB2, *ConfigCompUnbound().cmIntfDefinition) == false) {
      return false;
   }
   if (intU.read(pU1, cU1, cU2, *ConfigCompUnbound().cmIntfDefinition) == false) {
      return false;
   }
   DistUqentryCmp intCmp(debug);
   intCmp.setIgnoreSizeDiff(true);
   float intDist = intCmp(intB, intU);
   return intDist < u1MinIntDist;
}

/*----------------------------------------------------------------------------*/

static uint32_t countGaps(const Chain* ic, const Residues &sir, float gapMaxCtoNdist) {
   uint32_t numGaps = 0;

   for (uint32_t i = 1; i < ic->resis().size(); i++) {
      const Residue *prev = ic->resis()[i-1];
      const Residue *curr = ic->resis()[i];

      if (prev->c() == NULL) {
         continue;
      }
      if (curr->n() == NULL) {
         continue;
      }

      /* is there a gap? */
      if (prev->c()->dist(*curr->n()) > gapMaxCtoNdist) {
         /* does any of the two residues belong to the interface? */
         bool interface = false;
         interface = interface || std::binary_search(sir.begin(), sir.end(), prev) == true;
         interface = interface || std::binary_search(sir.begin(), sir.end(), curr) == true;
         /* if so, count a gap */
         if (interface == true) {
            numGaps++;
         }
      }
   }
   return numGaps;
}

/*----------------------------------------------------------------------------*/

static uint32_t countGaps(const Chains &ics, const Residues &ir, float gapMaxCtoNdist) {
   uint32_t numGaps = 0;
   Residues sir = ir;
   std::sort(sir.begin(), sir.end());
   /* sum up for all chains */
   for (uint32_t i = 0; i < ics.size(); i++) {
      numGaps += countGaps(ics[i], sir, gapMaxCtoNdist);
   }
   return numGaps;
}

/*----------------------------------------------------------------------------*/

#define MIN_BBATOM_DIST 0.8f
#define MIN_BBATOM_DIST_SQ (MIN_BBATOM_DIST * MIN_BBATOM_DIST)

/*----------------------------------------------------------------------------*/

static uint32_t chainClash(const std::vector<Atom*> &a, const std::vector<Atom*> &b) {
   uint32_t clashes = 0;
   for (uint32_t i = 0; i < a.size(); i++) {
      for (uint32_t j = 0; j < b.size(); j++) {
         if (a[i]->squareDist(*b[j]) < MIN_BBATOM_DIST_SQ) {
            clashes++;
         }
      }
   }

   return clashes;
}

/*----------------------------------------------------------------------------*/
#define ALLOWED_INCOMPLETE_RATIO 0.25f
static bool checkBackbones(const Protein &p, DebugOut &dbg) {
   for (uint32_t i = 0; i < p.chains().size(); i++) {
      uint32_t numIncomplete = 0;
      for (uint32_t j = 0; j < p.chains()[i]->resis().size(); j++) {
         if (p.chains()[i]->resis()[j]->hasCompleteBackbone() == false) {
            numIncomplete++;
         }
         if ((float)(numIncomplete) / ALLOWED_INCOMPLETE_RATIO > p.chains()[i]->resis().size()) {
            dbg("%s chain %c has %u/%u residues with complete backbone", p.name().c_str(), p.chains()[i]->name(), numIncomplete, p.chains()[i]->resis().size());
            return false;
         }
      }
   }
   return true;
}

/*----------------------------------------------------------------------------*/
#define MAX_CLASHING_ATOMS 20
static bool checkChainClashes(const Protein &p, DebugOut &dbg) {
   std::vector<std::vector<Atom*> > bbchains;

   /* grab all backbone atoms */
   for (uint32_t i = 0; i < p.chains().size(); i++) {
      bbchains.push_back(std::vector<Atom*>());
      for (uint32_t j = 0; j < p.chains()[i]->atoms().size(); j++) {
         if (p.chains()[i]->atoms()[j]->isBackBone() == true) {
            bbchains.back().push_back(p.chains()[i]->atoms()[j]);
         }
      }
   }

   /* 20 backbone atoms with pairwise distance < 0.8A */
   for (uint32_t i = 0; i < bbchains.size(); i++) {
      for (uint32_t j = i+1; j < bbchains.size(); j++) {
         uint32_t clashes = chainClash(bbchains[i], bbchains[j]);
         if (clashes > MAX_CLASHING_ATOMS) {
            dbg("%s: chain clash between %c and %c of size %u", p.name().c_str(), p.chains()[i]->name(), p.chains()[j]->name(), clashes);
            return false;
         }
      }
   }
   return true;
}

/*----------------------------------------------------------------------------*/

static uint32_t countClashes(const Residues& chainA, const Residues& chainB, float minCaDist) {
   uint32_t numOverlap = 0;
   for (uint32_t i = 0; i < chainA.size(); i++) {
      if (chainA[i]->hasCa() == false) {
         continue;
      }
      for (uint32_t j = 0; j < chainB.size(); j++) {
         if (chainB[j]->hasCa() == false) {
            continue;
         }
         if (chainB[j]->ca()->dist(*chainA[i]->ca()) < minCaDist) {
            /* overlap found - use next Ca of ChainA */
            numOverlap++;
         }
      }
   }

   return numOverlap;
}

/*----------------------------------------------------------------------------*/

static uint32_t countClashes(const Chains &cU1, const Chains &cB2, float minCaDist) {
   uint32_t numOverlap = 0;
   /* sum up for all chains */
   for (uint32_t i = 0; i < cU1.size(); i++) {
      for (uint32_t j = 0; j < cB2.size(); j++) {
         numOverlap += countClashes(cU1[i]->resis(), cB2[j]->resis(), minCaDist);
      }
   }
   return numOverlap;
}

static const char* trimFront(const char* s) {
   while (s[0] == ' ') {
      s++;
   }
   return s;
}
static std::string resToStr(const Residue* r) {
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

/*----------------------------------------------------------------------------*/

static void pdbInfo(const Protein &p, DebugOut &dbg) {
   dbg("Protein %s    chains:%lu resis:%lu atoms:%lu    nonresis:%lu",
         p.name().c_str(), p.chains().size(), p.resis().size(),
         p.atoms().size(), p.nonResis().size());
   for (uint32_t i = 0; i < p.chains().size(); i++) {
      dbg("  chain %c residues %lu %lu", p.chains()[i]->name(),
            p.chains()[i]->resis().size(), p.chains()[i]->atoms().size());
   }
}

/*----------------------------------------------------------------------------*/

static UnboundMatchResult alignUnbound(const std::string &id, const std::string &fnpBC, char nB1, char nB2, const std::string &fnpU1, char nU1) {
   UnboundMatchResult res;

   res.status = UnboundMatchResult::ERR_OK;

   /* read pdb files */
   Protein pBC;
   Protein pU1;
   {
      if (pBC.readPdb(fnpBC.c_str()) == false) {
         res.status = UnboundMatchResult::ERR_PDBFILE;
         res.errInfo = common::s_printf("error reading PDB file %s", fnpBC.c_str());
         return res;
      }
      if (pU1.readPdb(fnpU1.c_str()) == false) {
         res.status = UnboundMatchResult::ERR_PDBFILE;
         res.errInfo = common::s_printf("error reading PDB file %s", fnpU1.c_str());
         return res;
      }
   }


   /**
    * PBC  protein complex - protein 1 bound to protein 2
    * CB1  chain of protein1 bound to protein 2
    * CB2  chain of protein2 bound to protein 1
    * PU1  unbound protein 1
    * CU1  unbound chain of protein 1 sequence identical to CB1
    */


   /* select chains by their name */
   std::string cnB1(1, nB1);
   std::string cnB2(1, nB2);
   std::string cnU1(1, nU1);
   Chains cB1 = IChains::getChainsByName(pBC, cnB1);
   Chains cB2 = IChains::getChainsByName(pBC, cnB2);
   Chains cU1 = IChains::getChainsByName(pU1, cnU1);



   res.cChainsB1i = cnB1;
   res.cChainsB2i = cnB2;
   res.cChainsB1 = cnB1;
   res.cChainsB2x = cnB2;
   res.uChainsU1s = cnU1;
   res.uChainsU1x = cnU1;


   ChainConnectionInfo conInfo;
   SeqIdInfo seqInfo;
   ChainOverlapInfo ovlInfo(ConfigCompUnbound().u1ToB2CaClashDist);

   res.seedPBC = common::removeExtension(common::removePath(fnpBC));
   res.seedB1 = cnB1;
   res.seedB2 = cnB2;
   res.seedPU1 = common::removeExtension(common::removePath(fnpU1));
   res.seedU1 = cnU1;

   bool debug = true;
   /* to keep track of current thread */
   DebugOutNull _dbgNull;
   DebugOutPrefix _dbgPfx(common::s_printf("%s %s %s %s %s %s | ", id.c_str(), res.seedPBC.c_str(), res.seedB1.c_str(), res.seedB2.c_str(), res.seedPU1.c_str(), res.seedU1.c_str()));
   DebugOut *_dgb = NULL;
   if (cmd.var("debug").getBool() == true) {
      _dgb = &_dbgPfx;
   } else {
      _dgb = &_dbgNull;
      debug = false;
   }
   DebugOut &dbg = *_dgb;

   dbg.pushPrefix("PDBinfo");
   pdbInfo(pBC, dbg);
   pdbInfo(pU1, dbg);
   dbg.popPrefix();

   // TODO:
   // short chains filter should be used here
   // xtalcompunbound $PDBDATADIR/2amf1.pdb E K $PDBDATADIR/2ahr2.pdb A

   if (cB1.size() == 0) {
      res.status = UnboundMatchResult::ERR_CHAINSB1;
      res.errInfo = common::s_printf("ERR_CHAINSB1");
      return res;
   }
   if (cB2.size() == 0) {
      res.status = UnboundMatchResult::ERR_CHAINSB2;
      res.errInfo = common::s_printf("ERR_CHAINSB2");
      return res;
   }
   if (cU1.size() == 0) {
      res.status = UnboundMatchResult::ERR_CHAINSU1;
      res.errInfo = common::s_printf("ERR_CHAINSU1");
      return res;
   }
   cB1.resize(1);
   cB2.resize(1);
   cU1.resize(1);

   res.seedPBC = common::removeExtension(common::removePath(fnpBC));
   res.seedPU1 = common::removeExtension(common::removePath(fnpU1));




   /* check PDB structures for CA only chains */
   if (checkBackbones(pBC, dbg) == false) {
      res.status |= UnboundMatchResult::ERR_BACKBONEB;
      return res;
   }
   if (checkBackbones(pU1, dbg) == false) {
      res.status |= UnboundMatchResult::ERR_BACKBONEU1;
      return res;
   }



   /* check PDB structures for self clashes (1dfj has wrong biounit) */
   if (checkChainClashes(pBC, dbg) == false) {
      res.status |= UnboundMatchResult::ERR_CHAINCLASHES;
      return res;
   }
   if (checkChainClashes(pU1, dbg) == false) {
      res.status |= UnboundMatchResult::ERR_CHAINCLASHES;
      return res;
   }




   /* check interface size of B1:B2 */
   {
      uint32_t intsize = conInfo.getInterfaceSize(cB1.at(0), cB2.at(0),
            *ConfigCompUnbound().cmIntfDefinition);
      if (intsize < ConfigCompUnbound().minIntSize) {
         res.errInfo = common::s_printf("intsize %u", intsize);
         res.status |= UnboundMatchResult::ERR_NO_B1B2INT;
         return res;
      }
   }


   /* match chains from U1 to BC (providing a seed)*/
   dbg.pushPrefix("ChainMatch");
   ChainMatch cm(conInfo, seqInfo, pBC, cB1.at(0), cB2.at(0), pU1, cU1.at(0), dbg);
   dbg.popPrefix();

   /* check match */
   if (cm.good() == false) {
      res.status |= UnboundMatchResult::ERR_CM_BAD;
      return res;
   }


   ChainMatchExpansion cmx;
   dbg.pushPrefix("ChainMatchExpansionB2");
   cmx.expandB2(conInfo, seqInfo, pBC, pU1, cB2.at(0), cm, dbg);
   dbg.popPrefix();


   /* compute interface (multiple chains) */
   dbg.pushPrefix("MultiChainInterface");
   MultiChainInterface mint(cmx.B1(), cmx.B2(), conInfo, ConfigCompUnbound().minIntSize, *ConfigCompUnbound().cmIntfDefinition, dbg);
   dbg.popPrefix();


   assert(mint.cR().size() > 0 && "B1 needs to have an interface");
   assert(mint.cL().size() > 0 && "B2 needs to have an interface");

   /* check if DNA is within B interface region */
   if (CofactorDetector::checkDna(pBC, mint.cR(), mint.cL(), *ConfigCompUnbound().cmIntfDefinition, ConfigCompUnbound().minDnaDist) == false) {
      res.status |= UnboundMatchResult::ERR_DNA;
      return res;
   }



   /* superimposition and interface RMSD */
   Residues intB1; /* aligned interface */
   Residues intU1; /* aligned interface */
   {
      /* use a bigger interface to superimposition */
      ChainConnectionInfo ci;
      MultiChainInterface mint10(cmx.B1(), cmx.B2(), ci, ConfigCompUnbound().minIntSize, *ConfigCompUnbound().superimoseIntfDefinition, _dbgNull);


      /* set number of interface Ca */
      Residues caResis1;
      Residues caResis2;
      Residue::getCaResidues(mint10.resR(), caResis1);
      Residue::getCaResidues(mint10.resL(), caResis2);
      res.cNumI1Ca = caResis1.size();
      res.cNumI2Ca = caResis2.size();

#if ANALYZE_BFACTORS
      /* analyze interface bfactors */
      if (debug == true) {
         for (uint32_t i = 0; i < caResis1.size(); i++) {
            dbg("bfactor: %6.2f", caResis1[i]->ca()->tempFactor());
         }
         for (uint32_t i = 0; i < caResis2.size(); i++) {
            dbg("bfactor: %6.2f", caResis2[i]->ca()->tempFactor());
         }
      }
#endif

      /* compute interface alignment */
      res.uAlignedIntRatio = getAlignedInterface(cm.B1(), cm.U1(), mint10.intR(), intB1,
            intU1, dbg, res);
      if (res.uAlignedIntRatio < ConfigCompUnbound().minB1U1alignedInterfaceRatio) {
         res.errInfo = common::s_printf("ratio %6.3f of %6.3f", res.uAlignedIntRatio, ConfigCompUnbound().minB1U1alignedInterfaceRatio);
         res.status |= UnboundMatchResult::ERR_INTALGN;
         return res;
      }

      /* superimpose */
      std::vector<Atom*> cafrom;
      std::vector<Atom*> cato;
      getCaAtoms(intB1, cato, intU1, cafrom);
      KabschWrapper::superimpose(pU1.atoms(), cafrom, cato, &res.rotTrans);

      /* check interface RMSD */
      res.uIRmsd = KabschWrapper::getRmsd(cato, cafrom);
      if (res.uIRmsd > ConfigCompUnbound().maxU1IRmsd) {
         res.errInfo = common::s_printf("irmsd %6.3f of %6.3f", res.uIRmsd, ConfigCompUnbound().maxU1IRmsd);
         res.status |= UnboundMatchResult::ERR_IMSD;
         if (debug == true) {
            return res;
         }
      }
   }

   // TODO: pre-read and share lists
   std::vector<std::string> cofIgnoreList;
   const char* fn_cof_grp = getenv(ENV_COF_IGNORELIST);
   const char* fn_cof_ign = getenv(ENV_COF_GROUPS);
   cofIgnoreList = common::readList(fn_cof_ign);
   if (cofIgnoreList.empty()) {
      fprintf(stderr, "no entries in cofactor ignore list - check env var %s (\"%s\")\n", ENV_COF_GROUPS, fn_cof_grp);
      exit(1);
   }
   CofactorType ct;
   if (!ct.readGroupDef(fn_cof_grp)) {
      fprintf(stderr, "no entries in cofactor group definition - check env var %s (\"%s\")\n", ENV_COF_GROUPS, fn_cof_grp);
      exit(1);
   }
   dbg.setPrefix("CofactorDetector B1:");
   CofactorDetector cdB1(pBC, mint.resR(), ConfigCompUnbound().minCofactorDistB, ConfigCompUnbound().minCofactorNeighborsB, ConfigCompUnbound().cofNeighborDist, cofIgnoreList, dbg);
   dbg.setPrefix("CofactorDetector B2:");
   CofactorDetector cdB2(pBC, mint.resL(), ConfigCompUnbound().minCofactorDistB, ConfigCompUnbound().minCofactorNeighborsB, ConfigCompUnbound().cofNeighborDist, cofIgnoreList, dbg);
   dbg.setPrefix("CofactorDetector U1:");
   CofactorDetector cdU1(pU1, intU1, ConfigCompUnbound().minCofactorDistU, ConfigCompUnbound().minCofactorNeighborsU, ConfigCompUnbound().cofNeighborDist, cofIgnoreList, dbg);

   /* do not allow any DNA within interface of U1 */
   if (cdU1.hasDna() == true) {
      res.status |= UnboundMatchResult::ERR_DNAU1;
      return res;
   }

   dbg.setPrefix("CofactorMatcher B1->U1");
   CofactorMatcher cfm(CofactorType(), cdB1.cofactors(), intB1, cdU1.cofactors(), intU1, ct, dbg);
   dbg.unsetPrefix();
   UniqueCofactorDescriptor cduB1B2(cdB1, cdB2, cdU1, cfm);


   dbg.pushPrefix("ChainMatchExpansionU1");
   cmx.expandU1(conInfo, seqInfo, pBC, pU1, cB2.at(0), cm, dbg);
   dbg.popPrefix();


   if (cmx.UC().size() > 0 && ConfigCompUnbound().fixSmallPeptideClashes) {
      dbg.pushPrefix("ChainMatchExpansionU1fixclash");
      cmx.fixUCligand(dbg);
      dbg.popPrefix();
   }
   /* any conflicts with expansion? */
   if (cmx.UC().size() > 0 && ConfigCompUnbound().fixOligomerConflicts == true) {
      /* try to expand without seq-similars and see if that works */
      dbg.pushPrefix("ChainMatchExpansionU1oligo");
      cmx.expandU1oligo(conInfo, seqInfo, pBC, pU1, cB2.at(0), cm, dbg);
      dbg.popPrefix();
   }


   /* any clashes of cm.U1() with B2? */
   if (hasChainOverlap(cm.B1(), cB2.at(0), cm.U1(), ovlInfo,
         ConfigCompUnbound().u1ToB2MaxConflictRatio, dbg) == true) {
      res.status |= UnboundMatchResult::ERR_U1B2CLASH;
      if (debug == true) {
         return res;
      }
   }

   dbg.pushPrefix("ChainMatchExpansion");
   cmx.print(dbg);
   dbg.popPrefix();


   /* any chain of unbound in conflict with b2? */
   if (cmx.UC().empty() == false) {
      res.status |= UnboundMatchResult::ERR_CM_UC;
      if (debug == true) {
         return res;
      }
   }

   /* is there any interface in unbound U1:^U1 - similar to B1:B2? */
   if (hasIntSimilarity(pBC, cm.B1(), cmx.B2(), pU1, cm.U1(), ConfigCompUnbound().u1MinIntDist) == true) {
      res.status |= UnboundMatchResult::ERR_U1INTSIM;
      if (debug == true) {
         return res;
      }
   }

   /* disulfide bonds */
   std::vector<std::pair<Residue*, Residue*> > s2bonds;
   {
      s2bonds = Chain::getS2Bonds(mint.cR(), mint.cL(), 2.3f);
      for (uint32_t i = 0; i < s2bonds.size(); i++) {
         dbg("S2 bond %s %s", resToStr(s2bonds[i].first).c_str(), resToStr(s2bonds[i].second).c_str());
      }
   }


   res.cChainsB1i = chains2Str_sorted(mint.cR()).c_str();
   res.cChainsB2i = chains2Str_sorted(mint.cL()).c_str();
   res.cChainsB1 = chains2Str(cm.B1());
   res.cChainsB2x = chains2Str(cmx.B2());
   res.cNumGaps = 0;
   res.cNumI1Gaps = countGaps(mint.cR(), mint.resR(), ConfigCompUnbound().gapMaxCtoNdist);
   res.cNumI2Gaps = countGaps(mint.cL(), mint.resL(), ConfigCompUnbound().gapMaxCtoNdist);
   res.cNumGaps = res.cNumI1Gaps + res.cNumI2Gaps;
   res.cNumInterfaceChains = res.cChainsB1i.size() + res.cChainsB2i.size();
   res.cNumNonInterfaceChains = (res.cChainsB1.size() + res.cChainsB2x.size()) - res.cNumInterfaceChains;
   res.cNumLigands = cduB1B2.numCofactors();
   res.cNumS2Bonds = s2bonds.size();

   res.uChainsU1s = chains2Str(cm.U1());
   res.uChainsU1x = chains2Str(cmx.U1());


   res.u1NumChains = res.uChainsU1x.size();
   res.u1NumAddChains = res.uChainsU1x.size() - res.uChainsU1s.size();
   res.u1NumGaps = countGaps(cmx.U1(), intU1, ConfigCompUnbound().gapMaxCtoNdist);
   res.uNumMatchedLigands = cfm.numMatched();
   res.uNumUnmatchedLigands = cfm.numUnmatched();
   res.uNumClashes = countClashes(cmx.U1(), mint.cL(), ConfigCompUnbound().u1ToB2CaClashDist);

   res.cofStr = cduB1B2.cofStr();

   dbg("%s: %s  ", pBC.name().c_str(), pBC.pdbInfo().title);
   dbg("%s: %s", pU1.name().c_str(), pU1.pdbInfo().title);
   dbg("interface Ca RMSD of %s:%s and %s:%s is %.2fA", pBC.name().c_str(),
         res.cChainsB1.c_str(), pU1.name().c_str(), res.uChainsU1s.c_str(),
         res.uIRmsd);

   std::string sState = "";
   if (res.status == UnboundMatchResult::ERR_OK) {
      sState += "ok";
   } else {
      sState += "irmsd";
   }

   /* visualize */
#ifdef HAS_VIEWER   
   #pragma omp critical
   if (Viewer::getInstance().isRunning() == true) {
      Viewer::getInstance().remAll();
      Model *mdl;
      // set rotation and zoom
      {
         Orientator ori;
         ori.add(cmx.B1());
         ori.add(cmx.U1());
         ori.orient();
      }
      // focus on interface
      {
         Vector cen;
         cen += Residue::caCentre(mint.resR());
         cen += Residue::caCentre(mint.resL());
         cen /= 2;
         Viewer::getInstance().setCentre(cen);
      }

      for (uint32_t i = 0; i < cmx.B1().size(); i++) {
         RgbColor col = RgbColor(0x4488FF);
         mdl = new GlSseModel(*(cmx.B1()[i]), col);
         Viewer::getInstance().addModel(mdl);
      }

      for (uint32_t i = 0; i < cmx.B2().size(); i++) {
         RgbColor col = RgbColor(0x4444FF);
         mdl = new GlSseModel(*(cmx.B2()[i]), col);
         Viewer::getInstance().addModel(mdl);
      }

      for (uint32_t i = 0; i < cmx.U1().size(); i++) {
         RgbColor col = RgbColor(0x44FF44);
         mdl = new GlSseModel(*(cmx.U1()[i]), col);
         Viewer::getInstance().addModel(mdl);
      }

      for (uint32_t i = 0; i < cmx.UX().size(); i++) {
         RgbColor col = RgbColor(0xCCFFCC);
         mdl = new GlSseModel(*(cmx.UX()[i]), col);
         Viewer::getInstance().addModel(mdl);
      }

      for (uint32_t i = 0; i < cmx.UC().size(); i++) {
         RgbColor col = RgbColor(0xFFFF44);
         mdl = new GlSseModel(*(cmx.UC()[i]), col);
         Viewer::getInstance().addModel(mdl);
      }

      for (uint32_t i = 0; i < cmx.US().size(); i++) {
         RgbColor col = RgbColor(0xFFFF88);
         mdl = new GlSseModel(*(cmx.US()[i]), col);
         Viewer::getInstance().addModel(mdl);
      }

      for (uint32_t i = 0; i < cmx.BX().size(); i++) {
         RgbColor col = RgbColor(0x8888FF);
         mdl = new GlSseModel(*(cmx.BX()[i]), col);
         Viewer::getInstance().addModel(mdl);
      }

      for (uint32_t i = 0; i < cmx.BS().size(); i++) {
         RgbColor col = RgbColor(0xCC88FF);
         mdl = new GlSseModel(*(cmx.BS()[i]), col);
         Viewer::getInstance().addModel(mdl);
      }

      RgbColor col;
      col = RgbColor::DEEPPINK;
      col.t() = 0.25;
      mdl = new SticksModel(pU1.nonResis(), AtmColorer(col, RgbColor::PURPLE));
      Viewer::getInstance().addModel(mdl);

      SDL_Delay(50);
      cmd.setVar("viewer.screenshotfile", common::s_printf("/tmp/xtal_%s_%s-%s_%s_%s_%s.tga", sState.c_str(), pBC.name().c_str(), res.cChainsB1.c_str(), res.cChainsB2x.c_str(), pU1.name().c_str(), res.uChainsU1x.c_str()));
      Viewer::getInstance().screenshot();
   }
#endif   
	return res;
}


/*----------------------------------------------------------------------------*/

static bool startAlignUnbound(const std::string &id, const std::string &fnpBC, char nB1, char nB2,
      const std::string &fnpU1, char nU1, std::string &result) {
   UnboundMatchResult ret = alignUnbound(id, fnpBC, nB1, nB2, fnpU1, nU1);
   result = ret.toTableRow();
   return ret.status == UnboundMatchResult::ERR_OK;
}

/*----------------------------------------------------------------------------*/
#include <fstream>
template<class T>
static bool readTable(const std::string &filename, uint32_t nCols,
      std::vector<std::vector<T> > &table, bool append = false) {
   std::ifstream infile;

   if (nCols == 0) {
      return false;
   }
   if (append == false) {
      table.clear();
   }

   infile.open(filename.c_str(), std::ifstream::in);

   if (infile.good() == false) {
      return false;
   }

   std::vector<T> row(nCols);
   while (infile.good()) {
      for (uint32_t i = 0; i < row.size(); ++i) {
         infile >> row[i];
      }
      if (infile.good() == false) {
         break;
      }

      table.push_back(row);
   }
   infile.close();
   return table.size() != 0;
}

/*----------------------------------------------------------------------------*/

void XtalCompUnbound::registerStuff() {
   cmd.addVar("debug", "false");
}

/*----------------------------------------------------------------------------*/

int XtalCompUnbound::start() {
   int ret = 1;
   if (cmd.getNumArgs() == 5 || cmd.getNumArgs() == 6) {
      const std::string fnpBC = cmd.getArgStr(0);
      const char nB1 =  cmd.getArgStr(1).at(0);
      const char nB2 =  cmd.getArgStr(2).at(0);
      const std::string fnpU1 = cmd.getArgStr(3);
      const char nU1 =  cmd.getArgStr(4).at(0);

      std::string res;
      ret = startAlignUnbound("0", fnpBC, nB1, nB2, fnpU1, nU1, res) ? 0 : 1;
      printf("%s", res.c_str());
#ifdef HAS_VIEWER         
      /* auto mode? */
      bool autoMode = (cmd.getNumArgs() >= 6 && cmd.getArgStr(5) == "auto");
      if (Viewer::getInstance().isRunning() && autoMode == true) {
         Viewer::getInstance().stop(0);
      }
#endif      
   } else if (cmd.getNumArgs() == 2 || cmd.getNumArgs() == 3) {
      const std::string dnPdb = cmd.getArgStr(0);
      const std::string fnList = cmd.getArgStr(1);
      const std::string fnCheckpoint = cmd.getNumArgs() > 2 ? cmd.getArgStr(2) : std::string();

      std::vector<std::vector<std::string> > candList;
      readTable(fnList, 6, candList);
      /* remove header */
      candList.erase(candList.begin());

      if (candList.empty() == false) {
         ret = 0;
      }

#ifdef HAS_VIEWER     
      bool quit = false;
#endif      

      /* when viewing do not process parallel */
      //const uint32_t chunkSize = cmd.var("viewer").getBool() == true ? 1 : 1024;
      {
         ChunkStatus cs;
         cs.init(getName(), fnCheckpoint.c_str(), candList.size());
         cs.load();
         if (cs.getProgress() > 0) {
            Log::inf("resuming at %.3f%%", cs.getProgress() * 100);
         }
         #pragma omp parallel for
         for (uint32_t i = 0; i < candList.size(); i++) {
            if (cs.isCompleted(i)) {
               continue;
            }
            #pragma omp task
            {
               const std::string fnpBC = dnPdb + candList[i][1] + ".pdb";
               const char nB1 = candList[i][2].at(0);
               const char nB2 = candList[i][3].at(0);
               const std::string fnpU1 = dnPdb + candList[i][4] + ".pdb";
               const char nU1 = candList[i][5].at(0);
               std::string res;
               startAlignUnbound(candList[i][0].c_str(), fnpBC, nB1, nB2, fnpU1, nU1, res);
               #pragma omp critical
               {
                  printf("%s %s\n", candList[i][0].c_str(), res.c_str());
                  fflush(stdout);
                  cs.setCompleted(i);
                  if (i%1000 == 0) {
                     Log::inf("processed %6.3f%%", cs.getProgress() * 100);
                  }
               }
            }
         }
#ifdef HAS_VIEWER            
         if (cmd.var("viewer").getBool() == true && Viewer::getInstance().isRunning() == false) {
            quit = true;
         }
#endif         
      }

   } else {
      Log::err("usage: %s  PDB_BC CHAIN_B1 CHAIN_B2 PDB_U1 CHAIN_U1 [auto]", getName());
      Log::err("   checks aligns unbound structure to complex", getName());
#ifdef HAS_VIEWER         
      Viewer::getInstance().stop(0);
#endif
   }

   return ret;
}
