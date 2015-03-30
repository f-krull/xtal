#ifndef DEBUGOUTPUT_H_
#define DEBUGOUTPUT_H_

#include <stdarg.h>
#include "../libxtalutil/log.h"
#include <string>
#include <vector>
#include <stdint.h>

/*----------------------------------------------------------------------------*/

class DebugOut {
public:
   virtual ~DebugOut();
   virtual void operator ()(const char* msg, ...);

   virtual void setPrefix(const std::string &p);
   virtual void unsetPrefix();
   virtual void pushPrefix(const std::string &p);
   virtual void popPrefix();

protected:
   std::vector<std::string> m_prefixStack;
   std::string m_prefix;
   virtual void updatePrefix();
};

/*----------------------------------------------------------------------------*/

class DebugOutNull: public DebugOut {
public:
   virtual ~DebugOutNull() {};
   virtual void operator ()(const char* msg, ...);
private:
};

/*----------------------------------------------------------------------------*/

class DebugOutPrefix: public DebugOut {
public:
   DebugOutPrefix(const std::string &pfx);
   virtual ~DebugOutPrefix();
   virtual void operator ()(const char* msg, ...);
private:
   char *m_buffer;
   uint32_t m_n;
};

#endif /* DEBUGOUTPUT_H_ */

