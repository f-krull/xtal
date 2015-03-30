#ifndef INTFSEQCOMPARER_H_
#define INTFSEQCOMPARER_H_

#include <string>
#include <stdint.h>

/*----------------------------------------------------------------------------*/

class IntfSeqComparer {
public:
   IntfSeqComparer();

   /* interface is supposed to be represented as capitol letters */
   void align(const std::string &seq1, const std::string &seq2, uint32_t totalSeq, uint32_t totalInt);

   void alignMin(const std::string &seq1, const std::string &seq2);

   void label(const std::string &name1, const std::string &name2, char chain1 ,char chain2);

   float getIntScore() const;
   float getSeqId() const;
   float getSeqScore() const;
   void logAlignment(const void *p) const;

   void debug(const void *p) const;

   /* combine two similarities (a1,a2)(b1,b2) and (a1,a2)(b2,b1)
    * favor balanced deviation; seqsim <- function(x,y){(1-(1-x)^2) * (1-(1-y)^2)} */
   static float combineIntSim(const float seqsim1, const float seqsim2);

private:
   float m_interfaceScore;
   float m_seqScore;
   float m_seqId;

   uint32_t m_totalSeq;
   uint32_t m_totalInt;

   std::pair<std::string, std::string> m_alignment;

   std::string m_name1;
   std::string m_name2;
   char m_chain1;
   char m_chain2;
};

#endif /* INTFSEQCOMPARER_H_ */
