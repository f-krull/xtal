#include "iatoms.h"
#include <inttypes.h>
#include <errno.h>
#include <typeinfo>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include "../libxtalutil/log.h"

/*----------------------------------------------------------------------------*/

IAtoms::IAtoms() {
}

/*----------------------------------------------------------------------------*/

static inline void deleteAtom(Atom *a) {
   delete a;
}

/*----------------------------------------------------------------------------*/

#define READ_FIELD_STR(dst,src,offset,len) \
   if (offset + sizeof(dst)-1 < len) { \
      memcpy(dst, &src[offset], sizeof(dst)-1); \
      dst[sizeof(dst)-1] = '\0'; \
   } else { \
      memset(dst, 0, sizeof(dst)); \
   }

bool IAtoms::readPdb(const char* filename, PdbInfo* pdbinf) {
   FILE *file;
   char buffer[81];
   Atom a;

   enum STATE {
      PARSE_HEADER,
      PARSE_EXPDATA,
      PARSE_RESOLUTION,
      PARSE_ATOMS
   };

   STATE state = pdbinf != NULL ? PARSE_HEADER : PARSE_ATOMS;

   memset(pdbinf, 0, sizeof(PdbInfo));

   std::vector<Atom*> atoms;
   /* read atoms from files */
   file = fopen(filename, "r");
   if (file == NULL) {
      Log::err("%s: '%s' %s", typeid(*this).name(), filename, strerror(errno));
      return false;
   }
   while (fgets(buffer, sizeof(buffer) - 1, file) != NULL) {
#ifndef READ_MODELS
      /* stop at first model */
      if (strncmp(buffer, "ENDMDL", 6) == 0) {
         break;
      }
#endif
      if (state == PARSE_HEADER && strncmp(buffer, "HEADER", 6) == 0) {
         /* copy header data */
         READ_FIELD_STR(pdbinf->title, buffer, 10, sizeof(buffer)-1);
         READ_FIELD_STR(pdbinf->date, buffer, 50, sizeof(buffer)-1);
         READ_FIELD_STR(pdbinf->id, buffer, 62, sizeof(buffer)-1)
         /* read header - advance to next state */
         state = PARSE_EXPDATA;
      } else if (state == PARSE_EXPDATA && strncmp(buffer, "EXPDTA", 6) == 0) {
         READ_FIELD_STR(pdbinf->expdata, buffer, 10, sizeof(buffer)-1)
         state = PARSE_RESOLUTION;
      } else if (state == PARSE_RESOLUTION && strncmp(buffer, "REMARK   2 RESOLUTION.   ", 25) == 0) {
         READ_FIELD_STR(pdbinf->resolution, buffer, 25, sizeof(buffer)-1)
         state = PARSE_ATOMS;
      } else if (a.parse(buffer)) {
         atoms.push_back(new Atom(a));
         state = PARSE_ATOMS;
      }
   }
   fclose(file);
   if (atoms.size() == 0) {
      Log::err("%s: '%s' unable to read any atoms", typeid(*this).name(),
            filename);
   }
   IAtoms::destroy();
   IAtoms::set(atoms.begin(), atoms.end());
   return size() != 0;
}

/*----------------------------------------------------------------------------*/

void IAtoms::print() {
   for (size_t i = 0; i < atoms().size(); i++) {
      atoms()[i]->print();
   }
}

/*----------------------------------------------------------------------------*/

Atom *IAtoms::getByName(const char *name) const {
    for(uint32_t i = 0;i < atoms().size();i++){
        if(atoms()[i]->isName(name)){
            return atoms()[i];
        }
    }
    return NULL;
}

/*----------------------------------------------------------------------------*/

Vector IAtoms::centre() const {
   return Atom::centre(atoms());
}

/*----------------------------------------------------------------------------*/

void IAtoms::destroy() {
   std::for_each(m_atoms.begin(), m_atoms.end(), deleteAtom);
   m_atoms.clear();
   update();
}

/*----------------------------------------------------------------------------*/

IAtoms::~IAtoms() {

}

/*----------------------------------------------------------------------------*/

const std::vector<Atom*> & IAtoms::atoms() const {
   return m_atoms;
}

/*----------------------------------------------------------------------------*/

std::vector<Atom*> & IAtoms::atoms() {
   return m_atoms;
}

/*----------------------------------------------------------------------------*/

size_t IAtoms::size() const {
   return m_atoms.size();
}

/*----------------------------------------------------------------------------*/

void IAtoms::add(Atom *a) {
   m_atoms.push_back(a);
   update();
}

/*----------------------------------------------------------------------------*/

void IAtoms::update() {
}


