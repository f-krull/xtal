#include "multichaininterface.h"
#include "../libxtalcommon/kabschwrapper.h"
#include "../libxtalutil/log.h"
#include "configcompunbound.h"
#include "debugoutput.h"
#include <stdint.h>
#include <assert.h>

/*----------------------------------------------------------------------------*/

void MultiChainInterface::getInterface(const Chains &r, const Chains &l,
         Chains &cR, Chains &cL, std::vector<std::vector<cbool> > &intr,
         std::vector<std::vector<cbool> > &intl, ChainConnectionInfo &ci,
         uint32_t minIntfSize, const IntfDef &id, DebugOut &dbg) {

   {
      /* report connected chains only */
      std::set<Chain*> scR;
      std::set<Chain*> scL;

      dbg("chain interface sizes:");
      for (uint32_t i = 0; i < r.size(); i++) {
         for (uint32_t j = 0; j < l.size(); j++) {
            bool added = false;
            const uint32_t intsize = ci.getInterfaceSize(r[i], l[j], id);
            if (intsize >= minIntfSize) {
               scR.insert(r[i]);
               scL.insert(l[j]);
               added = true;
            }
            dbg("  %c - %c : %u  %s", r[i]->name(), l[j]->name(), intsize, added == true ? "(added)" : "");
         }
      }
      cR.assign(scR.begin(), scR.end());
      cL.assign(scL.begin(), scL.end());
   }

   /* init interface info */
   for (uint32_t i = 0; i < r.size(); i++) {
      intr.push_back(std::vector<cbool>(r[i]->resis().size(), false));
   }
   for (uint32_t j = 0; j < l.size(); j++) {
      intl.push_back(std::vector<cbool>(l[j]->resis().size(), false));
   }

   /* compare all chains */
   for (uint32_t i = 0; i < r.size(); i++) {
      /* new interface info r */
      for (uint32_t j = 0; j < l.size(); j++) {
         /* new interface info l */

         uint32_t cChainCont = 0;
         /* for every pair of chains compare all residues */
         for (uint32_t ii = 0; ii < r[i]->resis().size(); ii++) {
            for (uint32_t jj = 0; jj < l[j]->resis().size(); jj++) {
               /* residues already in contact? */
               if ((intr[i][ii] == true) && (intl[j][jj] == true)) {
                  continue;
               }
               if (id(r[i]->resis()[ii], l[j]->resis()[jj]) == true) {
                  intr[i][ii] = true;
                  intl[j][jj] = true;
                  cChainCont++;
               }
            }
         }
      }
   }

}

/*----------------------------------------------------------------------------*/

static void addResis(const std::vector<Chain*> &chains,
         const std::vector<std::vector<cbool> > &intf,
         std::set<Residue*> &resis) {
   assert(chains.size() == intf.size());
   for (uint32_t i = 0; i < chains.size(); ++i) {
      assert(chains[i]->resis().size() == intf[i].size());
      for (uint32_t j = 0; j < chains[i]->resis().size(); ++j) {
         if (intf[i][j] == true) {
            resis.insert(chains[i]->resis()[j]);
         }
      }
   }
}

/*----------------------------------------------------------------------------*/

MultiChainInterface::MultiChainInterface(const Chains &r, const Chains &l,
         ChainConnectionInfo &ci, uint32_t minIntfSize, const IntfDef &id, DebugOut &dbg) {
   getInterface(r, l, m_cR, m_cL, m_intr, m_intl, ci, minIntfSize, id, dbg);

   std::set<Residue*> resr;
   std::set<Residue*> resl;

   /* collect interface residue pointers */
   addResis(r, m_intr, resr);
   addResis(l, m_intl, resl);

   /* extract R Ca */
   m_ica.reserve(resr.size() + resl.size());
   for (std::set<Residue*>::iterator it = resr.begin(); it != resr.end();
            ++it) {
      if ((*it)->ca() != NULL) {
         m_ica.push_back((*it)->ca());
      }
   }
   for (std::set<Residue*>::iterator it = resl.begin(); it != resl.end();
            ++it) {
      if ((*it)->ca() != NULL) {
         m_ica.push_back((*it)->ca());
      }
   }
   /* make copies used as reference for RMSD */
   m_ica_cpy.resize(m_ica.size());
   for (uint32_t i = 0; i < m_ica.size(); i++) {
      m_ica_cpy[i] = new Atom(*m_ica[i]);
   }

   m_resr.assign(resr.begin(), resr.end());
   m_resl.assign(resl.begin(), resl.end());
}

/*----------------------------------------------------------------------------*/

MultiChainInterface::~MultiChainInterface() {
   for (uint32_t i = 0; i < m_ica_cpy.size(); ++i) {
      delete m_ica_cpy[i];
   }
}

/*----------------------------------------------------------------------------*/

std::vector<std::vector<cbool> > MultiChainInterface::intR() const {
   return m_intr;
}

/*----------------------------------------------------------------------------*/

std::vector<std::vector<cbool> > MultiChainInterface::intL() const {
   return m_intl;
}

/*----------------------------------------------------------------------------*/

const Residues& MultiChainInterface::resR() const {
   return m_resr;
}

/*----------------------------------------------------------------------------*/

const Residues& MultiChainInterface::resL() const {
   return m_resl;
}

/*----------------------------------------------------------------------------*/

const Chains& MultiChainInterface::cR() const {
   return m_cR;
}

/*----------------------------------------------------------------------------*/

const Chains& MultiChainInterface::cL() const {
   return m_cL;
}

/*----------------------------------------------------------------------------*/

float MultiChainInterface::getInterfaceRmsd() const {
   return KabschWrapper::getRmsd(m_ica, m_ica_cpy);
}
