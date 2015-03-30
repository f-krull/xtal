
#ifndef XTALUNIQUECOMP_H_
#define XTALUNIQUECOMP_H_

#include "../libxtalcommon/exectemplate.h"
#include <string>

/*----------------------------------------------------------------------------*/

class XtalUniqueCompPriv;
class UqEntry;

/*----------------------------------------------------------------------------*/

class XtalUniqueComp : public ExecTemplate {
public:
   XtalUniqueComp();
   virtual ~XtalUniqueComp();
protected:
private:
   XtalUniqueCompPriv *m;


   bool read(const std::string &filename, UqEntry &entry);
   const char* getName() const;
   int start();

};

#endif /*XTALUNIQUECOMP_H_*/
