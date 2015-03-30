#include "ichains.h"
#include <algorithm>
#include <assert.h>

/*----------------------------------------------------------------------------*/

IChains::IChains() {
}

/*----------------------------------------------------------------------------*/

//#include <stdio.h>

bool IChains::build(IAtoms &a, IResidues &r) {
   char a_lastchain;

   /* algorithm needs at least one atom to work */
   if (a.size() == 0) {
      return false;
   }

   /* initialize with first atom */
   std::vector<Atom*>::iterator ita_start = a.atoms().begin();
   /* initialize with first residue */
   std::vector<Residue*>::iterator itr_start = r.resis().begin();

   a_lastchain = (*ita_start)->chainId();

   std::vector<Atom*>::iterator ita_curr;
   std::vector<Residue*>::iterator itr_curr;

   ita_curr = ita_start;
   /* while there are atoms left */
   while (ita_curr != a.atoms().end()) {
      /* goto next atom */
      ita_curr++;

      /* if chain is new, search residue */
      if (ita_curr == a.atoms().end()
            || a_lastchain != (*ita_curr)->chainId()) {

         //printf("a] found chain %c\n", a_lastchain);

         /* search belonging residue */
         itr_curr = itr_start;
         while (itr_curr != r.resis().end()
               && (*itr_curr)->atoms()[0]->chainId() == a_lastchain) {
            itr_curr++;
         }

         chains().push_back(
               new Chain(ita_start, ita_curr, itr_start, itr_curr));
         if (itr_start != r.resis().end()
               && (*itr_start)->atoms()[0]->chainId() == a_lastchain) {
            //printf("b] found chain %c\n", a_lastchain);
         }

         /* set seed for continued search */
         itr_start = itr_curr;
         if (ita_curr != a.atoms().end()) {
            a_lastchain = (*ita_curr)->chainId();
         }
         ita_start = ita_curr;
      }

   }
   return true;
}

/*----------------------------------------------------------------------------*/

static inline void deleteChain(Chain *c) {
   delete c;
}

/*----------------------------------------------------------------------------*/

Vector IChains::centre() const {
   return Chain::centre(chains());
}

/*----------------------------------------------------------------------------*/

void IChains::destroy() {
   std::for_each(m_chains.begin(), m_chains.end(), deleteChain);
   m_chains.clear();
   update();
}

/*----------------------------------------------------------------------------*/

IChains::~IChains() {
}

/*----------------------------------------------------------------------------*/

const std::vector<Chain*> & IChains::chains() const {
   return m_chains;
}

/*----------------------------------------------------------------------------*/

std::vector<Chain*> & IChains::chains() {
   return m_chains;
}

/*----------------------------------------------------------------------------*/

size_t IChains::size() const {
   return m_chains.size();
}

/*----------------------------------------------------------------------------*/

void IChains::add(Chain *c) {
   m_chains.push_back(c);
   update();
}

/*----------------------------------------------------------------------------*/

void IChains::update() {
}

/*----------------------------------------------------------------------------*/

void IChains::selectByName(const std::string &chNames) const {
   m_selection.clear();
   Chain::getByName(m_chains, chNames, m_selection);
}

/*----------------------------------------------------------------------------*/

Chains IChains::getChainsByName(const IChains &p, const std::string &cNames) {
   std::vector<Chain*> selected;
   p.selectByName(cNames);
   selected.insert(selected.begin(), p.getSelection().begin(),
            p.getSelection().end());
   return selected;
}

/*----------------------------------------------------------------------------*/

const Chains& IChains::getSelection() const {
   return m_selection;
}
