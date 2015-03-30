

#ifndef COFACTORGROUPING_H_
#define COFACTORGROUPING_H_


#include "../libxtaldata/protein.h"
#include "debugoutput.h"
#include "cofactordetector.h"
#include <set>
/*----------------------------------------------------------------------------*/

class CofactorType {
public:
   enum MatchType {
     MATCH_NONE = 0,
     MATCH_BY_NAME = 1,
     MATCH_BY_GROUP = 2,
     MATCH_BY_ATOMS = 3
   };

   static const char* getMatchTypeStr(MatchType t);

   bool readGroupDef(const char* fn);

   MatchType matches(const Residue* c1, const Residue* c2) const;

private:

   bool matchByName(const Residue* c1, const Residue* c2) const;
   bool matchByGroup(const Residue* c1, const Residue* c2) const;
   bool matchByAtoms(const Residue* c1, const Residue* c2) const;

   std::vector<std::set<std::string> > m_groups;
};

/*----------------------------------------------------------------------------*/

typedef std::string CofType;

class CofFamily {
public:

   CofFamily(Residue* r);

   bool matches(const Residue* r, const CofactorType &c) const;

   void add(Residue* r);

   std::string type() const;

   std::string toStr() const;

   const Residues & cofactors() const;

private:

   std::set<CofType> m_types;
   Residues m_cofactors;
};

/*----------------------------------------------------------------------------*/

class CofGrouping {
public:

   CofGrouping(const CofactorType &cm, const Residues &cof);

   std::string toStr() const;

   const std::vector<CofFamily> & families() const;

private:
   std::vector<CofFamily> m_families;

};

/*----------------------------------------------------------------------------*/

#include <stdint.h>
#include <map>

class CofactorMatcher {
public:
   CofactorMatcher(const CofactorType &c, const std::vector<CofactorDetector::CofNeighborhood> &cofB, const Residues &intB, const std::vector<CofactorDetector::CofNeighborhood> &cofU, const Residues &intU, const CofactorType &ct, DebugOut &dbg);
   uint32_t numMatched() const;
   uint32_t numUnmatched() const;

   Residue* getMatched(Residue *res) const;

private:
   uint32_t m_numUnmatched;
   uint32_t m_numMatched;

   std::vector<uint32_t> m_matchedCofU;
   std::map<Residue*, Residue*> m_cofMatch;
};

#endif /* COFACTORGROUPING_H_ */
