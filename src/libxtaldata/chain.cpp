
#include "chain.h"
#include <stdint.h>

/*----------------------------------------------------------------------------*/

Chain::Chain(std::vector<Atom*>::iterator & ab,
      std::vector<Atom*>::iterator & aend, std::vector<Residue*>::iterator & rb,
      std::vector<Residue*>::iterator & rend) {
   if (ab != aend) {
      name() = (*ab)->chainId();
      IAtoms::set(ab, aend);
   }
   if (rb != rend) {
      IResidues::set(rb, rend);
   }
}

/*----------------------------------------------------------------------------*/

Chain::~Chain() {
}

/*----------------------------------------------------------------------------*/

Vector Chain::centre(const Chains & c) {
   Vector centre;
   uint32_t n = 0;

   for (unsigned int i = 0; i < c.size(); i++) {
      Vector cen = centre = c[i]->centre();
      if (cen != Vector(0,0,0)) {
         centre = c[i]->centre() + centre;
         n++;
      }
   }
   return centre /= n;
}

/*----------------------------------------------------------------------------*/

Vector Chain::centre() const {
   return IResidues::centre();
}

/*----------------------------------------------------------------------------*/

char & Chain::name() {
   return m_name;
}

/*----------------------------------------------------------------------------*/

char Chain::name() const {
   return m_name;
}

/*----------------------------------------------------------------------------*/

void Chain::getByName(const Chains &chains, const std::string &selection, Chains &selected) {
   selected.clear();
   for (uint32_t i = 0; i < chains.size(); i++) {
      for (uint32_t j = 0; j < selection.size(); j++) {
         if (chains[i]->name() == selection[j]) {
            selected.push_back(chains[i]);
         }
      }
   }
}

/*----------------------------------------------------------------------------*/

static bool hasS2Contact(const Residue *r1, const Residue *r2, float cutoff) {
   for (uint32_t i = 0; i < r1->atoms().size(); i++) {
      if (r1->atoms()[i]->isName(Atom::ATM_Sg) == false) {
         continue;
      }
      for (uint32_t j = 0; j < r2->atoms().size(); j++) {
         if (r2->atoms()[j]->isName(Atom::ATM_Sg) == false) {
            continue;
         }
         return r1->atoms()[i]->dist(*r2->atoms()[j]) < cutoff;
      }
   }
   return false;
}

/*----------------------------------------------------------------------------*/

std::vector<std::pair<Residue*, Residue*> > Chain::getS2Bonds(const Residues &c1, const Residues &c2, float cutoff) {
   std::vector<std::pair<Residue*, Residue*> > s2bonds;
   for (uint32_t i = 0; i < c1.size(); i++) {
      if (c1[i]->atoms().front()->isResn(Atom::RES_CYS) == false) {
         continue;
      }
      for (uint32_t j = 0; j < c2.size(); j++) {
         if (c2[j]->atoms().front()->isResn(Atom::RES_CYS) == false) {
            continue;
         }
         if (hasS2Contact(c1[i], c2[j], cutoff) == true) {
            s2bonds.push_back(std::make_pair(c1[i], c2[j]));
         }
      }
   }
   return s2bonds;
}

/*----------------------------------------------------------------------------*/

std::vector<std::pair<Residue*, Residue*> > Chain::getS2Bonds(const Residues &c1, float cutoff) {
   std::vector<std::pair<Residue*, Residue*> > s2bonds;
   for (uint32_t i = 0; i < c1.size(); i++) {
      if (c1[i]->atoms().front()->isResn(Atom::RES_CYS) == false) {
         continue;
      }
      for (uint32_t j = i+1; j < c1.size(); j++) {
         if (c1[j]->atoms().front()->isResn(Atom::RES_CYS) == false) {
            continue;
         }
         if (hasS2Contact(c1[i], c1[j], cutoff) == true) {
            s2bonds.push_back(std::make_pair(c1[i], c1[j]));
         }
      }
   }
   return s2bonds;
}


/*----------------------------------------------------------------------------*/

std::vector<std::pair<Residue*, Residue*> > Chain::getS2Bonds(const Chains &c1, const Chains &c2, float cutoff) {
   std::vector<std::pair<Residue*, Residue*> > s2bonds;
   for (uint32_t i = 0; i < c1.size(); i++) {
      for (uint32_t j = 0; j < c2.size(); j++) {
         std::vector<std::pair<Residue*, Residue*> > tmp = Chain::getS2Bonds(c1[i]->resis(), c2[j]->resis(), cutoff);
         s2bonds.insert(s2bonds.end(), tmp.begin(), tmp.end());
      }
   }
   return s2bonds;
}
