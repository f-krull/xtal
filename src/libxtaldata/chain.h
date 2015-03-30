#ifndef CHAIN_H_
#define CHAIN_H_

#include "iatoms.h"
#include "iresidues.h"

/*----------------------------------------------------------------------------*/

class Chain;

/*----------------------------------------------------------------------------*/

typedef std::vector<Chain*> Chains;

/*----------------------------------------------------------------------------*/

class Chain: public IAtoms, public IResidues {
public:
   Chain();
   Chain(const Chain &r2);
   Chain(Atoms::iterator &ab, Atoms::iterator &aend, Residues::iterator &rb,
         Residues::iterator &rend);
   virtual ~Chain();

   Vector centre() const;

   char name() const;

   char & name();

   static Vector centre(const Chains &c);

   static void getByName(const Chains &chains, const std::string &chNames, Chains &selected);

   static std::vector<std::pair<Residue*, Residue*> > getS2Bonds(const Residues &c1, float cutoff = 2.3f);
   static std::vector<std::pair<Residue*, Residue*> > getS2Bonds(const Residues &c1, const Residues &c2, float cutoff = 2.3f);
   static std::vector<std::pair<Residue*, Residue*> >  getS2Bonds(const Chains &c1, const Chains &c2, float cutoff= 2.3f);

private:
   char m_name;
};

#endif /* CHAIN_H_ */
