
#ifndef XTALCOMPSEQID_H_
#define XTALCOMPSEQID_H_

#include "../libxtalcommon/exectemplate.h"

class XtalCompSeqId : public ExecTemplate {
public:
   virtual ~XtalCompSeqId() {}
protected:
private:
   const char* getName() const;
   int start();
};

#endif /* XTALCOMPSEQID_H_ */
