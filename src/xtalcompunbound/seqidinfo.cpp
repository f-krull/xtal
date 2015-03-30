
#include "seqidinfo.h"
#include "../libxtalcommon/seqaligner.h"
#include <float.h>

/*----------------------------------------------------------------------------*/

float getSeqScore(const std::string &seq1, const std::string &seq2) {
   const UnitScore us = UnitScore();
   float score = SeqAligner::align(seq1.c_str(), seq2.c_str(), NULL, us, SeqAligner::SPEEDUP_NONE, '-');
   return score / us.getBestScore();
}

/*----------------------------------------------------------------------------*/

SeqIdInfo::SeqIdInfo() {
}

/*----------------------------------------------------------------------------*/

const SeqIdInfo::IdInf& SeqIdInfo::idInf(const Chain *b, const Chain *u, bool autoCompute) {
   std::map<std::pair<const Chain*, const Chain*>, IdInf>::const_iterator it;
   it = m_seqId.find(std::make_pair(b, u));
   /* not contained - compute */
   if (it == m_seqId.end()) {
      if (autoCompute == true) {
         const std::string sB = Residue::getResSequence(b->resis());
         const std::string sU = Residue::getResSequence(u->resis());
         if (sB.size() > 0 && sU.size() > 0) {
            float seqId = getSeqScore(sB, sU);
            m_sm.min = seqId / std::max(sB.size(), sU.size());
            m_sm.max = seqId / std::min(sB.size(), sU.size());
            m_sm.valid = true;
            //Log::dbg("SeqIdInfo: %c-%c %4.2f %f %4.2s", b->name(), u->name(), sm.min, sm.max, sm.valid ? "OK" : "");
         } else {
            m_sm.min = 0;
            m_sm.max = 0;
            m_sm.valid = false;
         }
         m_seqId[std::make_pair(b,u)] = m_sm;
         m_seqId[std::make_pair(u,b)] = m_sm;
         return m_seqId[std::make_pair(b,u)];
      } else {
         /* autoCompute == false */
         m_sm.min = 0;
         m_sm.max = 0;
         m_sm.valid = false;
         return m_sm;
      }
   }
   /* return interface size */
   return it->second;
}
