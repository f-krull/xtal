#ifndef MULTICHAININTERFACE_H_
#define MULTICHAININTERFACE_H_

#include "../libxtaldata/chain.h"
#include "debugoutput.h"
#include <set>
#include "chainconnectioninfo.h"

typedef char cbool;

/*----------------------------------------------------------------------------*/

class MultiChainInterface {
public:

   MultiChainInterface(const Chains &r, const Chains &l,
            ChainConnectionInfo &ci, uint32_t minIntfSize, const IntfDef &id, DebugOut &dbg);

   ~MultiChainInterface();

   std::vector<std::vector<cbool> > intR() const;
   std::vector<std::vector<cbool> > intL() const;

   /* we might want to know, which chains are involved in interface */
   const Chains& cR() const;
   const Chains& cL() const;

   const Residues& resR() const;
   const Residues& resL() const;

   float getInterfaceRmsd() const;

protected:

   static void getInterface(const Chains &r, const Chains &l, Chains &cR,
            Chains &cL, std::vector<std::vector<cbool> > &intR,
            std::vector<std::vector<cbool> > &intL, ChainConnectionInfo &ci,
            uint32_t minIntfSize, const IntfDef &id, DebugOut &dbg);

   Chains m_cR;
   Chains m_cL;

   std::vector<std::vector<cbool> > m_intr;
   std::vector<std::vector<cbool> > m_intl;

   std::vector<Residue*> m_resr;
   std::vector<Residue*> m_resl;

   std::vector<Atom*> m_ica;
   std::vector<Atom*> m_ica_cpy;
private:
};

#endif /* MULTICHAININTERFACE_H_ */
