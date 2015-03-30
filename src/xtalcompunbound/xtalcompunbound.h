
#ifndef XTALTEST_H_
#define XTALTEST_H_

#include "../libxtalcommon/exectemplate.h"

/*----------------------------------------------------------------------------*/

class XtalCompUnbound : public ExecTemplate {
public:
   virtual ~XtalCompUnbound() {}
   void registerStuff();
protected:
private:
   const char* getName() const;
   int start();
};

#endif /* XTALTEST_H_ */

