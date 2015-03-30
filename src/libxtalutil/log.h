#ifndef LOG_H_
#define LOG_H_

/*----------------------------------------------------------------------------*/

#include <stdarg.h>

class LoggerPriv;
class ILogObserver;

/*----------------------------------------------------------------------------*/

class Log {
public:
   static void dbg(const char *msg, va_list args);
   static void dbg(const char *msg, ...);
   static void war(const char *msg, ...);
   static void err(const char *msg, ...);
   static void inf(const char *msg, ...);

   static void addLogObserver(ILogObserver *o);
   static void remLogObserver(ILogObserver *o);

private:
   LoggerPriv *m;
   Log();
   ~Log();
   Log(const Log &);
   Log & operator=(const Log &);
   static Log& getInstance();
};

#endif /* LOG_H_ */
