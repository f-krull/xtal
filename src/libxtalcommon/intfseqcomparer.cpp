
#include "intfseqcomparer.h"
#include "../libxtalutil/log.h"
#include "../libxtalcommon/seqaligner.h"
#include <algorithm>
#include <assert.h>
#include <stdint.h>
#include <math.h>

#define XTALDEBUG

/*----------------------------------------------------------------------------*/

/* upper case and equal */
static bool isIntAligned(const char a, const char b) {
   if ((a == b) && isupper(a)) {
      return true;
   }
   return false;
}

/*----------------------------------------------------------------------------*/

class IntScore {
public:
   int32_t operator()(const char a, const char b) const {
      if ((a | 0x20) == (b | 0x20)) {
         return MAX_SCORE;
      }
      return -1;
   }

   int32_t max() const {
      return MAX_SCORE;
   }
private:
   static const int32_t MAX_SCORE;
};

const int32_t IntScore::MAX_SCORE = 2;


/*----------------------------------------------------------------------------*/

IntfSeqComparer::IntfSeqComparer() {
   m_interfaceScore = 0;
   m_seqScore = 0;
   m_seqId = 0;
   m_name1 = "";
   m_name2 = "";
   m_chain1 = ' ';
   m_chain2 = ' ';
   m_totalSeq = 0;
   m_totalInt = 0;
}

/*----------------------------------------------------------------------------*/

void IntfSeqComparer::alignMin(const std::string &seq1, const std::string &seq2) {
   uint32_t totalInt1 = std::count_if(seq1.begin(), seq1.end(), isupper);
   uint32_t totalInt2 = std::count_if(seq1.begin(), seq1.end(), isupper);
   uint32_t totalSeq1 = seq1.size();
   uint32_t totalSeq2 = seq1.size();
   align(seq1, seq2, std::min(totalSeq1, totalSeq2), std::min(totalInt1, totalInt2));
}

/*----------------------------------------------------------------------------*/

void IntfSeqComparer::align(const std::string &seq1, const std::string &seq2, uint32_t totalSeq, uint32_t totalInt) {
    m_totalSeq = totalSeq;
    m_totalInt = totalInt;


    /* do sequence alignment */
    m_seqScore = SeqAligner::align(seq1.c_str(), seq2.c_str(), &m_alignment, IntScore(), SeqAligner::SPEEDUP_SSE_32, '-');
    m_seqScore /= IntScore().max(); /* normalize to 1 score per res */

    /* count number aligned interface residues */
    uint32_t numAlignedInt = 0;
    assert(m_alignment.first.size() == m_alignment.second.size());
    for (uint32_t i = 0; i < m_alignment.first.size(); i++) {
       numAlignedInt += isIntAligned(m_alignment.first[i], m_alignment.second[i]) == true ? 1 : 0;
    }
    m_interfaceScore = numAlignedInt;
    m_interfaceScore = m_totalInt > 0 ? m_interfaceScore / m_totalInt : 0;

    assert(m_interfaceScore >= 0.);
    assert(m_interfaceScore <= 1.);

    if (m_totalSeq > 0) {
       m_seqId = ((float) m_seqScore / m_totalSeq);
    } else {
       m_seqId = 0;
    }
}

/*----------------------------------------------------------------------------*/

float IntfSeqComparer::getIntScore() const {
   return m_interfaceScore;
}

/*----------------------------------------------------------------------------*/

float IntfSeqComparer::getSeqScore() const {
   return m_seqScore;
}

/*----------------------------------------------------------------------------*/

float IntfSeqComparer::getSeqId() const {
   return m_seqId;
}

/*----------------------------------------------------------------------------*/

void IntfSeqComparer::logAlignment(const void *p) const {
   const char blank = ' ';
   Log::inf("%p alignment %s - %s", p, m_name1.c_str(), m_name2.c_str());
   Log::inf("%p seq1 %c    : %s", p, m_chain1, m_alignment.first.c_str());
   Log::inf("%p seq2 %c    : %s", p, m_chain2, m_alignment.second.c_str());
   std::string res;
   res.resize(m_alignment.first.size());
   assert(m_alignment.first.size() == m_alignment.second.size());
   for (uint32_t i = 0; i < res.size(); i++) {
      res[i] = IntScore()(m_alignment.first[i], m_alignment.second[i]) > 0 ? '*' : blank;
   }
   Log::inf("%p alignment : %s", p, res.c_str());
   for (uint32_t i = 0; i < res.size(); i++) {
      res[i] = isIntAligned(m_alignment.first[i], m_alignment.second[i]) == true ? '+' : blank;
   }
   Log::inf("%p interface : %s", p, res.c_str());
   for (uint32_t i = 0; i < res.size(); i++) {
      res[i] = isupper(m_alignment.first[i]) == 0 ? blank : 'i';
   }
   Log::inf("%p int1 %c    : %s", p, m_chain1, res.c_str());
   for (uint32_t i = 0; i < res.size(); i++) {
      res[i] = isupper(m_alignment.second[i]) == 0 ? blank : 'i';
   }
   Log::inf("%p int2 %c    : %s", p, m_chain2, res.c_str());
   Log::inf("%p %c->%c  int: %6.3f (%0.f/%u)  seq: %6.3f (%.0f/%u)", p, m_chain1, m_chain2, m_interfaceScore, m_interfaceScore * m_totalInt, m_totalInt, m_seqId, m_seqId * m_totalSeq, m_totalSeq);
}

/*----------------------------------------------------------------------------*/

void IntfSeqComparer::label(const std::string &name1, const std::string &name2, char chain1 ,char chain2) {
   m_name1 = name1;
   m_name2 = name2;
   m_chain1 = chain1;
   m_chain2 = chain2;
}

/*----------------------------------------------------------------------------*/

void IntfSeqComparer::debug(const void *p) const {
//   Log::inf("intscore: %6.3f  seqID: %6.3f (score: %.0f)  %s %s", getIntScore(), getSeqId(), getSeqScore(), m_name1.c_str(), m_name2.c_str());
//   if ((getSeqId() > 0.9 && getIntScore() < 0.6) || (getSeqId() < 0.6 && getIntScore() > 0.8)) {
    logAlignment(p);
//   }
}

/*----------------------------------------------------------------------------*/
/* favor balanced deviation; seqsim <- function(x,y){(1-(1-x)^2) * (1-(1-y)^2)} */
float IntfSeqComparer::combineIntSim(const float seqsim1, const float seqsim2) {
   return (1.f - pow(1.f - seqsim1, 2)) * (1.f - pow(1.f - seqsim2, 2));
}
