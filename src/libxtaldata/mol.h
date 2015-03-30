#ifndef MOL_H_
#define MOL_H_

#include "iatoms.h"
#include <vector>

class Mol: public IAtoms {
private:
protected:
public:

   Mol();
   virtual ~Mol();
   Mol(const Mol &m);

   bool readPdb(const char *filename);
   const std::string & name() const;
   std::string & name();
   const PdbInfo & pdbInfo() const;

private:
   std::string m_name;
   PdbInfo m_pdbInfo;
};

#endif /* MOL_H_ */
