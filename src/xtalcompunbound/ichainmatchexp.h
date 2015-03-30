#ifndef ICHAINMATCHEXP_H_
#define ICHAINMATCHEXP_H_

/*----------------------------------------------------------------------------*/

#include "../libxtaldata/chain.h"
#include "seqidinfo.h"
#include "configcompunbound.h"
#include "debugoutput.h"
#include <set>

/*----------------------------------------------------------------------------*/

class IChainMatchExp {
public:

   virtual ~IChainMatchExp()  {};

   const Chains & B1() const; /* matched pB1 chains */
   const Chains & U1() const; /* matched U1 or connected to matched U1 */
   const Chains & B2() const; /* B2 chain or pB2 chain directly or indirectly connected to B2 */
   /* unpleasant */
   const Chains & BX() const; /* not connected to B2 */
   const Chains & BS() const; /* connected to B2 but similar to some U1 */
   const Chains & UX() const; /* not connected to U1 or  */
   const Chains & UC() const; /* connected to U1 but in conflict with B2 */
   const Chains & US() const; /* connected to U1 but seq-similar to U1 */

   void print(DebugOut &dbg) const;

protected:

   Chains m_B1;
   Chains m_U1;
   Chains m_B2;

   Chains m_UX;
   Chains m_UC;
   Chains m_US;
   Chains m_BX;
   Chains m_BS;

   bool hasNoSimilar(const Chain *cx, const std::set<Chain*> &sCmU1,
            SeqIdInfo &si, float maxSeqID, DebugOut &dbg);

private:
};


#endif /* ICHAINMATCHEXP_H_ */
