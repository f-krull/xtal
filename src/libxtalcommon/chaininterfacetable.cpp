#include "chaininterfacetable.h"
#include "interface.h"
#include "../libxtalutil/log.h"
#include <assert.h>

/*----------------------------------------------------------------------------*/

ChainInterface::ChainInterface() {
   c1 = NULL;
   c2 = NULL;
   size1 = 0;
   size2 = 0;
   numBreaks = 0;
}

/*----------------------------------------------------------------------------*/

static uint32_t countInterfaceRes(const std::vector<cbool> &intf) {
   uint32_t numinterface;

   numinterface = 0;
   for (uint32_t i = 0; i < intf.size(); i++) {
      if (intf[i] == true) {
         numinterface++;
      }
   }
   return numinterface;
}

/*----------------------------------------------------------------------------*/

ChainInterface::ChainInterface(const Chain *c1, const Chain *c2, const IntfDef &id) {
   this->c1 = c1;
   this->c2 = c2;
   /* compute interface residues */
   Interface::getInterface(c1->resis(), c2->resis(), int1, int2, &id);
   /* count size of interfaces */
   size1 = countInterfaceRes(int1);
   size2 = countInterfaceRes(int2);
   numBreaks = 0;
}

/*----------------------------------------------------------------------------*/

void ChainInterface::flip() {
   const Chain *c = c1;
   c1 = c2;
   c2 = c;

   int1.swap(int2);
   /* numBreaks = numBreaks */
   uint32_t sizeTmp = size1;
   size1 = size2;
   size2 = sizeTmp;
}

/*----------------------------------------------------------------------------*/

ChainInterfaceTable::ChainInterfaceTable(const IntfDef &id) :
         m_interfaceCutoff(id) {
}

/*----------------------------------------------------------------------------*/

ChainInterfaceTable::~ChainInterfaceTable() {
   std::set<ChainInterface*>::iterator it = m_infUniq.begin();
   while (it != m_infUniq.end()) {
      delete *it;
      it++;
   }
}

/*----------------------------------------------------------------------------*/

static inline std::pair<const Chain*, const Chain*> getChainKey(const Chain* c1,
         const Chain *c2) {
   return std::make_pair(c1, c2);
}

/*----------------------------------------------------------------------------*/

ChainInterface * ChainInterfaceTable::add(const Chain *c1, const Chain *c2) {
   ChainInterface *intf_f;
   ChainInterface *intf_r;
   std::map<std::pair<const Chain*, const Chain*> , ChainInterface*>::const_iterator it;
   it = m_intfTable.find(getChainKey(c1, c2));
   if (it != m_intfTable.end()) {
      return it->second;
   }
   intf_f = new ChainInterface(c1, c2, m_interfaceCutoff);
   m_intfTable[getChainKey(c1, c2)] = intf_f;
   /* reverse order for reverse key*/
   intf_r = new ChainInterface();
   *intf_r = *intf_f;
   intf_r->flip();
   m_intfTable[getChainKey(c2, c1)] = intf_r;
   /* to be able to delete later on */
   m_infUniq.insert(intf_f);
   m_infUniq.insert(intf_r);
   return intf_f;
}

/*----------------------------------------------------------------------------*/

void ChainInterfaceTable::add(const Chains &sel, const Chains &all) {
   Chain *chain1, *chain2;
   for (uint32_t i = 0; i < sel.size(); i++) {
      for (uint32_t j = 0; j < all.size(); j++) {
         chain1 = sel[i];
         chain2 = all[j];
         /* avoid self interface */
         if (chain1 == chain2) {
            continue;
         }
         add(chain1, chain2);
      }
   }
}

/*----------------------------------------------------------------------------*/

ChainInterface * ChainInterfaceTable::get(const Chain *c1, const Chain *c2) {
   ChainInterface* ret = getCached(c1,c2);
   if (ret == NULL) {
      ret = add(c1, c2);
   }
   return ret;
}

/*----------------------------------------------------------------------------*/

ChainInterface* ChainInterfaceTable::getCached(const Chain *c1, const Chain *c2) {
   std::map<std::pair<const Chain*, const Chain*> , ChainInterface*>::iterator it;

   it = m_intfTable.find(getChainKey(c1, c2));
   if (it != m_intfTable.end()) {
      return it->second;
   }
   return NULL;
}

/*----------------------------------------------------------------------------*/

const ChainInterface* ChainInterfaceTable::getCached(const Chain *c1, const Chain *c2) const {
   std::map<std::pair<const Chain*, const Chain*> , ChainInterface*>::const_iterator it;

   it = m_intfTable.find(getChainKey(c1, c2));
   if (it != m_intfTable.end()) {
      return it->second;
   }
   return NULL;
}


/*----------------------------------------------------------------------------*/
/* mark additional interface residues */
static void addInt(std::vector<cbool> &to, const std::vector<cbool> &from) {
   assert(to.size() == from.size());
   for (uint32_t i = 0; i < to.size(); i++) {
      to.at(i) = to.at(i) || from.at(i);
   }
}

/*----------------------------------------------------------------------------*/
/* highlight all residues with interface from chain c to any of chain other */
std::string ChainInterfaceTable::getSeqInt(const Chain* c, const std::vector<const Chain*> &other) {
   std::string cseq = Residue::getResSequence(c->resis());
   std::vector<cbool> cintf(cseq.size(), false);

   for (uint32_t i = 0; i < other.size(); i++) {
      const ChainInterface *intf = this->get(c, other[i]);
      assert(intf != NULL);
      addInt(cintf, intf->int1);
   }
   for (uint32_t i = 0; i < cintf.size(); i++) {
      cseq[i] = cintf.at(i) ? cseq[i] : tolower(cseq[i]);
   }

   return cseq;
}

/*----------------------------------------------------------------------------*/

std::string ChainInterfaceTable::getSeqInt(const Chain* c, const std::vector<Chain*> &other) {
   std::vector<const Chain*> o(other.begin(), other.end());
   return getSeqInt(c, o);
}
