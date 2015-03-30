#ifndef CHAININTERFACETABLE_H_
#define CHAININTERFACETABLE_H_

#include "../libxtaldata/chain.h"
#include "../libxtalcommon/interface.h"
#include <stdint.h>
#include <set>
#include <map>

typedef char cbool;

/*----------------------------------------------------------------------------*/

class ChainInterface {
public:
   ChainInterface(const Chain *c1, const Chain *c2, const IntfDef &id);
   ChainInterface();
   const Chain *c1;
   const Chain *c2;
   uint32_t size1;
   uint32_t size2;
   std::vector<cbool> int1;
   std::vector<cbool> int2;
   uint32_t numBreaks;

   void flip();
private:
};

/*----------------------------------------------------------------------------*/

#if 0
/* cache sequence of chain */
class ChainSeqInfo {
public:

   const std::string& getSeq(const Chain *c) {
      std::map<const Chain*, std::string>::const_iterator it = m_seqs.find(c);
      if (it != m_seqs.end()) {
         return it->second;
      }
      m_seqs[c] = Residue::getResSequence(c->resis());
      return m_seqs[c];
   }

private:
   std::map<const Chain*, std::string> m_seqs;
};
#endif

/*----------------------------------------------------------------------------*/

class ChainInterfaceTable {
public:
   ChainInterfaceTable(const IntfDef &id);
   ~ChainInterfaceTable();

   void add(const Chains &sel, const Chains &all);
   ChainInterface * add(const Chain *c1, const Chain *c2);

   ChainInterface *get(const Chain *c1, const Chain *c2);
   ChainInterface *getCached(const Chain *c1, const Chain *c2);
   const ChainInterface *getCached(const Chain *c1, const Chain *c2) const;


   std::string getSeqInt(const Chain* c, const std::vector<Chain*> &other);
   std::string getSeqInt(const Chain* c, const std::vector<const Chain*> &other);

private:
   /* keeps track of any interfaces to delete them later on */
   std::set<ChainInterface*> m_infUniq;
   /* lookup table for all interfaces; two chains -> interface */
   std::map<std::pair<const Chain*, const Chain*>, ChainInterface*> m_intfTable;
   const IntfDef &m_interfaceCutoff;
};

#endif /* CHAININTERFACETABLE_H_ */
