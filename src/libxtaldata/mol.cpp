#include "mol.h"
#include "../libxtalutil/log.h"
#include "../libxtalutil/common.h"
#include <stdio.h>
#include <errno.h>
#include <stdint.h>
#include <string.h>
#include <typeinfo>
#include <algorithm>

/*----------------------------------------------------------------------------*/

Mol::Mol() {
}

/*----------------------------------------------------------------------------*/

Mol::~Mol() {
   IAtoms::destroy();
}

/*----------------------------------------------------------------------------*/

Mol::Mol(const Mol &m) {
   IAtoms::destroy();
   std::vector<Atom*> cpy(m.atoms().size());
   for (uint32_t i = 0; i < m.atoms().size(); i++) {
      cpy[i] = new Atom(*m.atoms()[i]);
   }
   IAtoms::set(cpy.begin(), cpy.end());
   m_name = m.m_name;
}

/*----------------------------------------------------------------------------*/

const std::string & Mol::name() const {
    return m_name;
}

/*----------------------------------------------------------------------------*/

std::string & Mol::name() {
    return m_name;
}

/*----------------------------------------------------------------------------*/

const PdbInfo & Mol::pdbInfo() const {
    return m_pdbInfo;
}

/*----------------------------------------------------------------------------*/

bool Mol::readPdb(const char *filename) {
   bool ret = IAtoms::readPdb(filename, &m_pdbInfo);
   /* assign PDB ID from filename if no header was found during reading */
   if (ret && (isdigit(m_pdbInfo.id[0]) == 0) && filename != NULL) {
      m_name = common::removePath(filename);
      m_name = common::removeExtension(m_name);
   } else {
      m_name = m_pdbInfo.id;
   }
   return ret;
}
