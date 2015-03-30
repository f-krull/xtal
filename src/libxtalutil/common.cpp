#include "common.h"
#include <fstream>
#include <stdarg.h>

/*----------------------------------------------------------------------------*/

void common::wait(unsigned int ms) {
   //SDL_Delay(ms);
}

/*----------------------------------------------------------------------------*/

bool common::exists(const std::string filename) {
   std::ifstream file;

   file.open(filename.c_str(), std::ios::in);
   if (file == NULL) {
      return false;
   }
   file.close();
   return true;
}

/*----------------------------------------------------------------------------*/

std::string common::s_printf(const char *str, ...) {
   char buffer[1024];

   va_list args;
   va_start(args, str);
   vsnprintf(buffer, sizeof(buffer)-1, str, args);
   va_end(args);
   buffer[sizeof(buffer)-1] = '\0';
   return std::string(buffer);
}

/*----------------------------------------------------------------------------*/

std::vector<std::string> common::readList(const char * filename, char comment) {
   std::vector<std::string> ret;
   std::ifstream file;
   std::string line;

   file.open(filename);
   if (file != NULL) {
      while (getline((file), line)) {
         if ((line.length() > 0) && (line.at(0) != comment)) {
            ret.push_back(line);
         }
      }
   }
   file.close();
   return ret;
}

/*----------------------------------------------------------------------------*/

std::string common::removeExtension(const std::string & filename) {
   std::string res = filename;
   // cut before last "."
   size_t pos_p = filename.find_last_of('.');
   if (pos_p != filename.npos) {
      res = res.substr(0, pos_p);
   }
   return res;
}

/*----------------------------------------------------------------------------*/

std::string common::removePath(const std::string & filename) {
   std::string res = filename;
   // cut at last "/"
   size_t pos_p = filename.find_last_of('/');
   if (pos_p != filename.npos && (pos_p + 1) < res.size()) {
      res = res.substr((pos_p + 1), res.size() - (pos_p + 1));
   }
   return res;
}

