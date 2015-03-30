#include "decoyreader.h"
#include "decoyreaderq.h"
#include "decoyreadermv.h"
#include "decoyreaderzd.h"
#include <string.h>

DecoyReader* DecoyReader::getDecoyReader(Mol *rec, Mol *lig, const char* filename) {
   const char *suffix;

   suffix = strrchr(filename, '.');
   if (suffix != NULL) {
      if (strcmp(suffix, ".mv") == 0) {
         return new DecoyReaderMv(rec, lig, filename);
      } else if (strcmp(suffix, ".dec") == 0) {
         return new DecoyReaderQ(rec, lig, filename);
      } else if (strcmp(suffix, ".out") == 0) {
         return new DecoyReaderZd(rec, lig, filename);
      }
   }
   return new DecoyReaderQ(rec, lig, filename);
}
