#ifndef KABSCHWRAPPER_H_
#define KABSCHWRAPPER_H_

#include "../libxtaldata/mol.h"
#include <vector>

/*----------------------------------------------------------------------------*/

class KabschWrapper {
protected:
public:
   static float getCaRmsd(const Mol &a, const Mol &b);
   static float getRmsd(const std::vector<Atom*> &a,
         const std::vector<Atom*> &b);

   static void superimpose(std::vector<Atom*> &atoms, const std::vector<
         Atom*> &from, const std::vector<Atom*> &to, std::vector<float> *mv = NULL);
};

#endif /*KABSCHWRAPPER_H_*/
