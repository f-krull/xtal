#ifndef ICHAINS_H_
#define ICHAINS_H_

#include "chain.h"

/*----------------------------------------------------------------------------*/

class IChains {
public:
   virtual ~IChains();

   const virtual Chains & chains() const;

   virtual Chains & chains();

   virtual size_t size() const;

   virtual bool build(IAtoms &atoms, IResidues &resis);

   Vector centre() const;

   virtual void selectByName(const std::string &chNames) const;

   virtual const Chains& getSelection() const;

   static Chains getChainsByName(const IChains &p, const std::string &cNames);


protected:
   IChains();

   template<typename T> void set(T begin, T end);

   template<typename T> void insert(T begin, T end);

   void add(Chain*);

   void destroy();

   virtual void update();

protected:
   Chains m_chains;
   mutable Chains m_selection;
};

#endif /* ICHAINS_H_ */
