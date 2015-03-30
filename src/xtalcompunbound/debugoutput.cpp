#include "debugoutput.h"
#include <stdio.h>
/*----------------------------------------------------------------------------*/

DebugOut::~DebugOut() {
}

/*----------------------------------------------------------------------------*/

void DebugOut::operator ()(const char* msg, ...) {
   va_list args;
   va_start(args, msg);
   Log::dbg(msg);
   va_end(args);
}

/*----------------------------------------------------------------------------*/

void DebugOut::setPrefix(const std::string &p) {
   m_prefix = p;
}

/*----------------------------------------------------------------------------*/

void DebugOut::unsetPrefix() {
   setPrefix("");
}

/*----------------------------------------------------------------------------*/

void DebugOut::updatePrefix() {
   m_prefix = "";
   for (uint32_t i = 0; i < m_prefixStack.size(); i++) {
      m_prefix += m_prefixStack[i];
   }
}

/*----------------------------------------------------------------------------*/

void DebugOut::pushPrefix(const std::string &p) {
   m_prefixStack.push_back(p);
   updatePrefix();
}

/*----------------------------------------------------------------------------*/

void DebugOut::popPrefix() {
   if (m_prefixStack.size() > 0) {
      m_prefixStack.pop_back();
      updatePrefix();
   }
}

/*----------------------------------------------------------------------------*/

void DebugOutNull::operator ()(const char* msg, ...) {
}

/*----------------------------------------------------------------------------*/
#define LOGBUFSIZE (1024 * 1024)
DebugOutPrefix::DebugOutPrefix(const std::string &pfx) {
   m_buffer = new char[LOGBUFSIZE];
   m_n = snprintf(m_buffer, LOGBUFSIZE - 1, "%s", pfx.c_str());
}

/*----------------------------------------------------------------------------*/

DebugOutPrefix::~DebugOutPrefix() {
   delete [] m_buffer;
}

/*----------------------------------------------------------------------------*/

void DebugOutPrefix::operator ()(const char* msg, ...) {
   uint32_t n2 = snprintf(m_buffer + m_n, LOGBUFSIZE - m_n, "%s%s", m_prefix.c_str(), m_prefix.empty() ? "" : " ");
   va_list args;
   va_start(args, msg);
   vsnprintf(m_buffer + m_n + n2, LOGBUFSIZE - m_n - n2, msg, args);
   va_end(args);
   Log::dbg(m_buffer);
}
