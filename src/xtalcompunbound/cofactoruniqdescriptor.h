#ifndef COFACTORUNIQDESCRIPTOR_H_
#define COFACTORUNIQDESCRIPTOR_H_

#include "cofactordetector.h"
#include "cofactorgrouping.h"

/*----------------------------------------------------------------------------*/

class UniqueCofactorDescriptor {
public:
   /* combines two detectors to generate one unique cofactor string for a complex */
   UniqueCofactorDescriptor(const CofactorDetector &dB1, const CofactorDetector &dB2, const CofactorDetector &dU1, const CofactorMatcher &match);

   const std::string & cofStr() const;

   uint32_t numCofactors() const;

private:

   std::string m_cofStr;
   uint32_t m_numCofactors;
};

#endif /* COFACTORUNIQDESCRIPTOR_H_ */
