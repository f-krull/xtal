#include "xtalcompseqid.h"
#include "../libxtalutil/log.h"
#include "../libxtalutil/command.h"
#include "../libxtalutil/common.h"
#include "../libxtaldata/protein.h"
#include "../libxtalcommon/seqaligner.h"
#include "../libxtalcommon/chaininterfacetable.h"

#include <stdio.h>
#include <assert.h>

#define CFG_SEQ_MIN_LENGTH 20
#define MIN_CHAIN_SIM 1.0f
#define SEQ_MIN_ID 0.250000001f

#define CFG_INTCUTOFF 6.0f
#define CFG_MINUMINTFRES 11

/*----------------------------------------------------------------------------*/

void filterShortChains(const IChains &chains, Chains &res) {
   res.clear();
   for (uint32_t i = 0; i < chains.size(); i++) {
      if (Residue::getResSequence(chains.chains()[i]->resis()).size() >= CFG_SEQ_MIN_LENGTH) {
         res.push_back(chains.chains()[i]);
      }
   }
}

/*----------------------------------------------------------------------------*/
/* if chain is unnamed, we can not discriminate it from others -> ignore */
void filterUnnamedChains(Chains &res) {
   Chains::iterator iter = res.begin();
   while (iter != res.end()) {
      /* space stand for unnamed */
      if ((*iter)->name() == ' ') {
         iter = res.erase(iter);
      } else {
         ++iter;
      }
   }
}

/*----------------------------------------------------------------------------*/

const char* XtalCompSeqId::getName() const {
	return "XtalCompSeqId";
}

/*----------------------------------------------------------------------------*/
class ChainSeqDescr {
public:
   std::string pdbId;
   std::vector<char> chainGroup;
   std::string seq;

   ChainSeqDescr(const std::string &pname, const Chain* c) {
      pdbId = pname;
      chainGroup.push_back(c->name());
      seq = Residue::getResSequence(c->resis());
   }

   char getGroupChain() const {
      if (chainGroup.size() > 0) {
         return chainGroup[0];
      }
      return '\0';
   }

   void mergeWith(const ChainSeqDescr &c) {
      chainGroup.insert(chainGroup.end(), c.chainGroup.begin(),
            c.chainGroup.end());
   }
};

/*----------------------------------------------------------------------------*/

float getSeqScore(const std::string &seq1, const std::string &seq2, const UnitScore &us) {
   float score = SeqAligner::align(seq1.c_str(), seq2.c_str(), NULL, us, SeqAligner::SPEEDUP_NONE, '-');
   return score / us.getBestScore();
}

/*----------------------------------------------------------------------------*/

uint32_t addChainDescr(const Protein &p, std::vector<ChainSeqDescr> &seqs, bool merge = true, bool debug = false) {
   std::vector<ChainSeqDescr> seqsTmp;

   Chains pchains;

   /* neglect short chains */
   filterShortChains(p, pchains);

   /* filter unnamed chains */
   filterUnnamedChains(pchains);

   for (uint32_t i = 0; i < pchains.size(); i++) {
      seqsTmp.push_back(ChainSeqDescr(p.name(), pchains[i]));
   }

   std::vector<ChainSeqDescr> seqsGroups;

   UnitScore us;

   /* group similar chains */
   bool added;
   for (uint32_t i = 0; i < seqsTmp.size(); i++) {
      added = false;

      if (merge == true) {
         /* try to match to previous clusters */
         for (uint32_t j = 0; j < seqsGroups.size(); j++) {
            /* compute the (worst) similarity */
            if (getSeqScore(seqsTmp[i].seq, seqsGroups[j].seq, us)
                  / std::max(seqsTmp[i].seq.size(), seqsGroups[j].seq.size()) >= MIN_CHAIN_SIM) {
               /* add similar chain to group */
               seqsGroups[j].mergeWith(seqsTmp[i]);
               added = true;
               break;
            }
         }
      }
      /* make new group */
      if (added == false) {
         seqsGroups.push_back(seqsTmp[i]);
      }
   }
   seqs.insert(seqs.end(), seqsGroups.begin(), seqsGroups.end());

   uint32_t maxLen = 0;
   for (uint32_t i = 0; i < seqsGroups.size(); i++) {
      if (seqsGroups[i].seq.size() > maxLen) {
         maxLen = seqsGroups[i].seq.size();
      }
   }
   if (debug == true) {
      Log::inf("%s chains %3lu valid %3lu merged %3lu (max len %u)", p.name().c_str(), p.chains().size(), seqsTmp.size(), seqsGroups.size(), maxLen);
   }
   return p.chains().size();
}

/*----------------------------------------------------------------------------*/

void compareSeq(const ChainSeqDescr &csd1, const ChainSeqDescr &csd2,
      const UnitScore &us, FILE *f, bool cutoff = true) {
   float score = getSeqScore(csd1.seq, csd2.seq, us);
   /* normalize by max score */

   /* normalize by smaller sequence length */
   float seqSim1 = score / csd1.seq.size();
   if (cutoff == false || seqSim1 >= SEQ_MIN_ID) {
      #pragma omp critical (sim_print)
      {
         fprintf(f, "%s %c %s %c %5.3f\n", csd1.pdbId.c_str(),
               csd1.chainGroup.at(0), csd2.pdbId.c_str(), csd2.chainGroup.at(0),
               seqSim1);
      }
   }
   float seqSim2 = score / csd2.seq.size();
   if (cutoff == false || seqSim2 >= SEQ_MIN_ID) {
      #pragma omp critical (sim_print)
      {
         fprintf(f, "%s %c %s %c %5.3f\n", csd2.pdbId.c_str(),
               csd2.chainGroup.at(0), csd1.pdbId.c_str(), csd1.chainGroup.at(0),
               seqSim2);
      }
   }
}

/*----------------------------------------------------------------------------*/

static bool getChainSeqDescr(const std::string &fnPdbPath,
      const std::vector<std::string> &pdbCodeList,
      std::vector<ChainSeqDescr> &csds, bool merge = false) {
   csds.clear();
   uint64_t numChains = 0;
   /* read pdb files */
   for (uint32_t i = 0; i < pdbCodeList.size(); i++) {
      Protein p;
      /* check if we were able to read */
      if (p.readPdb((fnPdbPath + "/" + pdbCodeList[i] + ".pdb").c_str())
            == false) {
         csds.clear();
         return false;
      }
      /* use filename as name (pdb5) */
      p.name() = pdbCodeList[i];
      numChains += addChainDescr(p, csds, merge);
   }
   Log::inf("successfully read %lu pdb files (%lu chains %lu valid sequences)",
         pdbCodeList.size(), numChains, csds.size());
   return true;
}

/*----------------------------------------------------------------------------*/

class PdbCodeLists {
public:
   FILE *fList1;
   FILE *fList2;
   FILE *fConf;
   const std::string m_fnOutPath;

   PdbCodeLists(std::string fnOutPath) :
         m_fnOutPath(fnOutPath) {
      fList1 = NULL;
      fList2 = NULL;
      fConf = NULL;
   }

   void close() {
      if (fList1 != NULL) {
         fclose(fList1);
      }
      if (fList2 != NULL) {
         fclose(fList2);
      }
      fList1 = NULL;
      fList2 = NULL;
   }

   bool open(uint32_t u) {
      std::string fnList1 = common::s_printf("c%04u.list1.txt", u);
      std::string fnList2 = common::s_printf("c%04u.list2.txt", u);
      /* open lists */
      fList1 = fopen((m_fnOutPath + "/" + fnList1).c_str(), "w");
      fList2 = fopen((m_fnOutPath + "/" + fnList2).c_str(), "w");
      /* open config file */
      fConf = fopen(
            common::s_printf("%s/c%04u.conf", m_fnOutPath.c_str(), u).c_str(),
            "w");

      if (fList1 == NULL || fList2 == NULL || fConf == NULL) {
         close();
         return false;
      }

      /* write config file */
      fprintf(fConf, "LIST1=%s\n", fnList1.c_str());
      fprintf(fConf, "LIST2=%s\n", fnList2.c_str());
      fclose(fConf);

      return true;
   }
};

/*----------------------------------------------------------------------------*/

static int writeConfList() {
   assert(cmd.getNumArgs() == 3);
   assert(cmd.getArgStr(0) == "split");

   const std::string fnInPdbList = cmd.getArgStr(1);
   const std::string fnOutPath = cmd.getArgStr(2);


   std::vector<std::string> pdbCodes = common::readList(fnInPdbList.c_str());

   if (pdbCodes.size() == 0) {
      Log::err("no pdbcodes could be read from: %s", fnInPdbList.c_str());
      return 1;
   }


   uint64_t countFiles = 0;
   const uint32_t entriesPerList1 = 1000;
   const uint32_t entriesPerList2 = 1000;
   PdbCodeLists lists(fnOutPath);

   for (uint32_t i = 0; i < pdbCodes.size(); i+=entriesPerList1) {
      for (uint32_t j = i+1; j < pdbCodes.size(); j+= entriesPerList2) {

         if (lists.open(countFiles) == false) {
            Log::err("unable to write to %s", fnOutPath.c_str());
            return 1;
         }
         countFiles++;
         assert(lists.fList1 != NULL);
         assert(lists.fList2 != NULL);

         uint32_t end1 = std::min((uint32_t) pdbCodes.size(),
               i + entriesPerList1);
         uint32_t end2 = std::min((uint32_t) pdbCodes.size(),
               j + entriesPerList2);


         for (uint32_t k = i; k < end1; k++) {
            fprintf(lists.fList1, "%s\n", pdbCodes[k].c_str());
         }

         for (uint32_t l = j; l < end2; l++) {
            fprintf(lists.fList2, "%s\n", pdbCodes[l].c_str());
         }

         lists.close();
      }
      Log::inf("wrote %6.3f%%", 100.0f * i / pdbCodes.size());
   }
   lists.close();
   return 0;
}

/*----------------------------------------------------------------------------*/

static bool cmpSim(const std::string &fnPdbPath,
      const std::vector<std::string> &list1,
      const std::vector<std::string> &list2, FILE *outfile, bool merge) {
   /* read PDB lists */
   std::vector<ChainSeqDescr> seqs1;
   std::vector<ChainSeqDescr> seqs2;

   bool ok = true;
   #pragma omp critical (sim_read1)
   {
      ok = ok && getChainSeqDescr(fnPdbPath, list1, seqs1, merge);
   }
   #pragma omp critical (sim_read2)
   {
      ok = ok && getChainSeqDescr(fnPdbPath, list2, seqs2, merge);
   }
   if (ok == false) {
      return false;
   }

   UnitScore us;

   uint64_t countTotal = 0;
   for (uint32_t i = 0; i < seqs1.size(); i++) {
      for (uint32_t j = 0; j < seqs2.size(); j++) {
         if (merge == true && (seqs1[i].pdbId < seqs2[j].pdbId) == false) {
            continue;
         }
         countTotal++;
      }
   }

   uint64_t countCurr = 0;
   for (uint32_t i = 0; i < seqs1.size(); i++) {
      for (uint32_t j = 0; j < seqs2.size(); j++) {
         if (merge == true && (seqs1[i].pdbId < seqs2[j].pdbId) == false) {
            continue;
         }
         compareSeq(seqs1[i], seqs2[j], us, outfile, merge);
         countCurr++;
//        if (countCurr % 1000 == 0) {
//           Log::inf("progress: %6.2f%% (%lu/%lu) chains",
//                 100.0f * countCurr / countTotal, countCurr, countTotal);
//        }
      }
   }
   return true;
}

/*----------------------------------------------------------------------------*/

static bool cmdChainCon(const std::string &fnPath, const std::string &fnList) {
   /* read paths to pdb files */
   std::vector<std::string> pdbCodeList = common::readList(fnList.c_str());

   /* make sure we read something */
   if (pdbCodeList.size() == 0) {
      Log::err("error reading: %s", fnList.c_str());
      return false;
   }

   #pragma omp parallel for
   for (uint32_t i = 0; i < pdbCodeList.size(); i++) {
      Protein p;
      /* check if we were able to read */
      if (p.readPdb((fnPath + "/" + pdbCodeList[i] + ".pdb").c_str()) == false) {
         continue;
      }

      /* use filename as name (pdb5) */
      p.name() = pdbCodeList[i];

      IntfDefAllAtom iaa(CFG_INTCUTOFF);
      ChainInterfaceTable citab(iaa);
      /* create con table */
      Chains pchains;
      filterShortChains(p, pchains);
      for (uint32_t k = 0; k < pchains.size(); k++) {
         for (uint32_t l = 0; l < pchains.size(); l++) {
            /* no chains interfacing themselves */
            if (k == l) {
               continue;
            }

            ChainInterface* cint = citab.add(pchains[k], pchains[l]);

            if (cint == NULL) {
               continue;
            }

            uint32_t intSize = std::max(cint->size1, cint->size2);

            if (intSize < CFG_MINUMINTFRES) {
               continue;
            }

            #pragma omp critical (sim_print)
            {
               printf("%s %c %c %u\n", p.name().c_str(), pchains[k]->name(),
                     pchains[l]->name(), intSize);
            }
         }
      }
   }
   return true;
}

/*----------------------------------------------------------------------------*/

static bool cmdChainGrp(const std::string &path, const std::string &list) {
   /* read paths to pdb files */
   std::vector<std::string> pdbCodeList = common::readList(list.c_str());

   /* make sure we read something */
   if (pdbCodeList.size() == 0) {
      Log::err("error reading: %s", list.c_str());
      return false;
   }

   #pragma omp parallel for
   for (uint32_t i = 0; i < pdbCodeList.size(); i++) {
      Protein p;
      /* check if we were able to read */
      if (p.readPdb((path + "/" + pdbCodeList[i] + ".pdb").c_str()) == false) {
         continue;
      }

      /* use filename as name (pdb5) */
      p.name() = pdbCodeList[i];

      /* create group table */
      std::vector<ChainSeqDescr> csds;
      addChainDescr(p, csds);
      for (uint32_t i = 0; i < csds.size(); i++) {
         for (uint32_t j = 0; j < csds[i].chainGroup.size(); j++) {
            #pragma omp critical (sim_print)
            {
               printf("%s %c %c\n", csds[i].pdbId.c_str(), csds[i].chainGroup[j],
                     csds[i].getGroupChain());
            }
         }
      }

   }
   return true;
}

/*----------------------------------------------------------------------------*/


static bool cmdChainSim(const std::string &path, const std::string &list) {
   /* read paths to pdb files */
   typedef std::vector<std::string> PdbCodes;

   PdbCodes pdbCodes = common::readList(list.c_str());

   /* make sure we read something */
   if (pdbCodes.size() == 0) {
      Log::err("error reading: %s", list.c_str());
      return false;
   }

   
   /* split to chunks of 2000x2000; use at least 4 chunks on smaller sets */
   const uint32_t entriesPerList = std::min(2000u, (uint32_t)(pdbCodes.size() / 4));

   std::vector<std::pair<PdbCodes, PdbCodes> > cmpLists;

   /* generate lists */
   for (uint32_t i = 0; i < pdbCodes.size(); i += entriesPerList) {
      for (uint32_t j = i + 1; j < pdbCodes.size(); j += entriesPerList) {
         cmpLists.push_back(
               std::make_pair(std::vector<std::string>(),
                     std::vector<std::string>()));

         uint32_t end1 = std::min((uint32_t) pdbCodes.size(),
               i + entriesPerList);
         uint32_t end2 = std::min((uint32_t) pdbCodes.size(),
               j + entriesPerList);

         for (uint32_t k = i; k < end1; k++) {
            cmpLists.back().first.push_back(pdbCodes[k]);
         }

         for (uint32_t l = j; l < end2; l++) {
            cmpLists.back().second.push_back(pdbCodes[l]);
         }

      }
      Log::inf("prepared %6.3f%%", 100.0f * i / pdbCodes.size());
   }

   /* process lists */
   uint32_t numDone = 0;
   bool ret = true; /* wont work yet with omp */
   #pragma omp parallel for
   for (uint32_t i = 0; i < cmpLists.size(); i++) {
      #pragma omp task shared(numDone)
      {
         ret = ret && cmpSim(path, cmpLists[i].first, cmpLists[i].second, stdout, true);
         #pragma omp atomic
         numDone++;
         Log::inf("processed %6.3f%%", 100.0f * numDone / cmpLists.size());
      }
   }
   return ret;
}

/*----------------------------------------------------------------------------*/

#define CMD_CON "con"
#define CMD_GRP "grp"
#define CMD_SIM "sim"

int XtalCompSeqId::start() {
   if (cmd.getNumArgs() < 3 || cmd.getNumArgs() > 5) {
		Log::err("usage: %s [PDBPATH] [PDB_LIST] [PDB_LIST] [OUTFILE]", getName());
		Log::err("   writes all pairwise chain sequence ids to OUTFILE");
      Log::err("usage: %s preanalyze [PDBPATH] [PDB_LIST] [OUTFILE1] [OUTFILE2]", getName());
		Log::err("   writes chain groups to OUTFILE1 and chain connectivity to OUTFILE2");
		Log::err("usage: %s split [PDB_LIST] [OUTDIR]", getName());
		Log::err("   splits a list of PDB codes to chunks");


		Log::err("usage: %s %s [PDBPATH] [PDB_LIST]", getName(), CMD_CON);
		Log::err("usage: %s %s [PDBPATH] [PDB_LIST]", getName(), CMD_GRP);
		Log::err("usage: %s %s [PDBPATH] [PDB_LIST]", getName(), CMD_SIM);

		return 1;
	}


   if (cmd.getArgStr(0) == CMD_CON) {
      const std::string path = cmd.getArgStr(1);
      const std::string list = cmd.getArgStr(2);
      return cmdChainCon(path, list) ? 0 : 1;
   }

   if (cmd.getArgStr(0) == CMD_GRP) {
      const std::string path = cmd.getArgStr(1);
      const std::string list = cmd.getArgStr(2);
      return cmdChainGrp(path, list) ? 0 : 1;
   }

   if (cmd.getArgStr(0) == CMD_SIM) {
      const std::string path = cmd.getArgStr(1);
      const std::string list = cmd.getArgStr(2);
      return cmdChainSim(path, list) ? 0 : 1;
   }

   if (cmd.getNumArgs() == 3) {
      return writeConfList();
   }

   const std::string fnPdbPath = cmd.getArgStr(0);
   const std::string fnPdbList1 = cmd.getArgStr(1);
   const std::string fnPdbList2 = cmd.getArgStr(2);
   const std::string fnOutfile = cmd.getArgStr(3);

   bool merge = fnOutfile != "clus_out.txt";

   /* open output file */
   FILE *outfile = fopen(fnOutfile.c_str(), "w");
   if (outfile == NULL) {
      Log::err("error opening output file %s", fnOutfile.c_str());
   }
   std::vector<std::string> list1 = common::readList(fnPdbList1.c_str());
   if (fnPdbList1.empty() == true) {
      Log::err("error reading list1 %s", fnPdbList1.c_str());
   }
   std::vector<std::string> list2 = common::readList(fnPdbList2.c_str());
   if (fnPdbList2.empty() == true) {
      Log::err("error reading list2 %s", fnPdbList2.c_str());
   }
   cmpSim(fnPdbPath, list1, list2, outfile, merge);
   fclose(outfile);
	return 0;
}
