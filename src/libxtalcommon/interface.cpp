#include "interface.h"
#include "kabschwrapper.h"
#include <stdint.h>

/*----------------------------------------------------------------------------*/

IntfDefAllAtom::IntfDefAllAtom(float cutoff) : m_cutoff(cutoff) {
}

/*----------------------------------------------------------------------------*/

cbool IntfDefAllAtom::operator()(const Residue *r, const Residue* o) const {
   return Residue::hasContact_fast(r, o, m_cutoff);
}

/*----------------------------------------------------------------------------*/

cbool IntfDefMendez2003::operator()(const Residue *r, const Residue* o) const {
   return Residue::hasContact_fast(r, o, 10.0f);
}

/*----------------------------------------------------------------------------*/

Interface::Interface(const std::vector<Residue*> &pr, const std::vector<
      Residue*> &pl, const IntfDef &id) {
   for (uint32_t i = 0; i < pr.size(); i++) {
      for (uint32_t j = 0; j < pl.size(); j++) {
         if (id(pr[i], pl[j]) == true) {
            m_ir.insert(pr[i]);
            m_il.insert(pl[j]);
         }
      }
   }
   /* extract R Ca */
   m_ica.reserve(m_ir.size() + m_il.size());
   for (std::set<Residue*>::iterator it = m_ir.begin(); it != m_ir.end(); ++it) {
      if ((*it)->ca() != NULL) {
         m_ica.push_back((*it)->ca());
      }
   }
   for (std::set<Residue*>::iterator it = m_il.begin(); it != m_il.end(); ++it) {
      if ((*it)->ca() != NULL) {
         m_ica.push_back((*it)->ca());
      }
   }
   /* make copies used as reference for RMSD */
   m_ica_cpy.resize(m_ica.size());
   for (uint32_t i = 0; i < m_ica.size(); i++) {
      m_ica_cpy[i] = new Atom(*m_ica[i]);
   }
}

/*----------------------------------------------------------------------------*/

Interface::~Interface() {
   for (uint32_t i = 0; i < m_ica_cpy.size(); ++i) {
      delete m_ica_cpy[i];
   }
}

/*----------------------------------------------------------------------------*/

std::vector<Residue*> Interface::getIntR() const {
   return std::vector<Residue*>(m_ir.begin(), m_ir.end());
}

/*----------------------------------------------------------------------------*/

std::vector<Residue*> Interface::getIntL() const {
   return std::vector<Residue*>(m_il.begin(), m_il.end());
}

/*----------------------------------------------------------------------------*/

float Interface::getInterfaceRmsd() const {
   return KabschWrapper::getRmsd(m_ica, m_ica_cpy);
}

/*----------------------------------------------------------------------------*/

void Interface::getInterface(const std::vector<Residue*> &chain1,
         const std::vector<Residue*> &chain2, std::vector<cbool> &int1,
         std::vector<cbool> &int2, const IntfDef *id) {
   int1.resize(chain1.size(), false);
   int2.resize(chain2.size(), false);
   for (uint32_t i = 0; i < chain1.size(); i++) {
      for (uint32_t j = 0; j < chain2.size(); j++) {
         if ((int1[i] == true) && (int2[j] == true)) {
            continue;
         }
         if ((*id)(chain1[i], chain2[j]) == true) {
            int1[i] = true;
            int2[j] = true;
         }
      }
   }
}

/*----------------------------------------------------------------------------*/

uint32_t Interface::getInterfaceSize(const Residues & r1, const Residues & r2,
         const IntfDef &id) {
   Interface intf(r1, r2, id);
   return intf.getIntR().size() * intf.getIntL().size();
}


