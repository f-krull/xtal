#include "log.h"
#include "ilogobserver.h"
#include <stdio.h>
#include <set>

/*----------------------------------------------------------------------------*/
#define LOGBUFSIZE (1024 * 1024)

#define LOG(msg) \
   va_list args; \
   va_start(args, msg); \
   vsnprintf(getInstance().m->buffer, LOGBUFSIZE - 1, msg, \
         args); \
   va_end(args); \
   getInstance().m->print();

/*----------------------------------------------------------------------------*/

class LoggerPriv {
private:
protected:
public:
   LoggerPriv() {
      buffer = new char[LOGBUFSIZE];
   }
   ~LoggerPriv() {
      delete [] buffer;
   }
   char *buffer;
   void print();
   std::set<ILogObserver*> obs;
};

/*----------------------------------------------------------------------------*/

void LoggerPriv::print() {
   fprintf(stderr, "%s\n", buffer);

   std::set<ILogObserver*>::iterator it;
   for (it = obs.begin(); it != obs.end(); it++) {
      (*it)->getLog(buffer);
   }
}

/*----------------------------------------------------------------------------*/

Log::Log() {
   m = new LoggerPriv();
}

/*----------------------------------------------------------------------------*/

Log::~Log() {
   delete m;
}

/*----------------------------------------------------------------------------*/

void Log::dbg(const char *msg, va_list args) {
   #pragma omp critical (log)
   {
      vsnprintf(getInstance().m->buffer, LOGBUFSIZE - 1, msg, args);
      getInstance().m->print();
   }
}

/*----------------------------------------------------------------------------*/

void Log::dbg(const char* msg, ...) {
   #pragma omp critical (log)
   {
      LOG(msg)
   }
}

/*----------------------------------------------------------------------------*/

void Log::inf(const char* msg, ...) {
   #pragma omp critical (log)
   {
      LOG(msg)
   }
}

/*----------------------------------------------------------------------------*/

void Log::war(const char* msg, ...) {
   #pragma omp critical (log)
   {
      LOG(msg)
   }
}

/*----------------------------------------------------------------------------*/

void Log::err(const char* msg, ...) {
   #pragma omp critical (log)
   {
      LOG(msg)
   }
}

/*----------------------------------------------------------------------------*/

Log& Log::getInstance() {
   static Log log;
   return log;
}

/*----------------------------------------------------------------------------*/

void Log::addLogObserver(ILogObserver *o) {
   getInstance().m->obs.insert(o);
}

/*----------------------------------------------------------------------------*/

void Log::remLogObserver(ILogObserver *o) {
   std::set<ILogObserver*>::iterator it;

   it = getInstance().m->obs.find(o);
   if (it != getInstance().m->obs.end()) {
      getInstance().m->obs.erase(it);
   }
}
