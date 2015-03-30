#include "residue.h"
#include "aamap.h"
#include "../libxtaldata/matrix.cpp"
#include "../libxtalutil/log.h"
#include "../libxtalutil/common.h"
#include <sstream>
#include <iomanip>
#include <climits>
#include <cstdlib>
#include <string.h>

/*----------------------------------------------------------------------------*/

Residue::Bond::Bond(Atom *a, Atom *b) {
   this->atm1 = a;
   this->atm2 = b;
}

/*----------------------------------------------------------------------------*/

char Residue::HELIX = 'H';
char Residue::SHEET = 'E';
char Residue::COIL = 'C';
char Residue::UNSET = '\0';

/*----------------------------------------------------------------------------*/

Residue::Residue() {
}

/*----------------------------------------------------------------------------*/

void Residue::init() {
   m_ssetype = UNSET;
   m_ca = NULL;
   m_cb = NULL;
   m_c = NULL;
   m_o = NULL;
   m_n = NULL;
}

/*----------------------------------------------------------------------------*/

Residue::~Residue() {
   for (unsigned int i = 0; i < m_bonds.size(); i++) {
      delete m_bonds[i];
   }
}

/*----------------------------------------------------------------------------*/

Residue::Residue(const Residue &r2) {
   atoms().resize(r2.atoms().size());
   for (unsigned int i = 0; i < r2.atoms().size(); i++) {
      atoms()[i] = new Atom((*r2.atoms()[i]));
   }
   this->m_ssetype = r2.m_ssetype;
   this->m_ca = r2.m_ca;
   this->m_cb = r2.m_cb;
   this->m_c = r2.m_c;
   this->m_o = r2.m_o;
   this->m_n = r2.m_n;
}

/*----------------------------------------------------------------------------*/

Residue::Residue(Atoms::iterator &begin, Atoms::iterator &end) {
   init();
   for ( ; begin != end; begin++) {
      assert((*begin) != NULL);
      /* check atom type and insert atom to res */
      if (strncmp((*begin)->name(), Atom::ATM_CA, 4) == 0) {
         m_ca = (*begin);
      } else if (strncmp((*begin)->name(), Atom::ATM_C, 4) == 0) {
         m_c = (*begin);
      } else if (strncmp((*begin)->name(), Atom::ATM_O, 4) == 0) {
         m_o = (*begin);
      } else if (strncmp((*begin)->name(), Atom::ATM_N, 4) == 0) {
         m_n = (*begin);
      } else if (strncmp((*begin)->name(), Atom::ATM_CB, 4) == 0) {
         m_cb = (*begin);
      }
      /* add to atoms */
      IAtoms::add((*begin));
   }
}

/*----------------------------------------------------------------------------*/

Residue & Residue::operator=(const Residue &r2) {
   atoms().resize(r2.atoms().size());
   for (unsigned int i = 0; i < r2.atoms().size(); i++) {
      atoms()[i] = new Atom((*r2.atoms()[i]));
   }
   this->m_ssetype = r2.m_ssetype;
   this->m_ca = r2.m_ca;
   this->m_cb = r2.m_cb;
   this->m_c = r2.m_c;
   this->m_o = r2.m_o;
   this->m_n = r2.m_n;
   return (*this);
}

/*----------------------------------------------------------------------------*/

const Atom* Residue::ca() const {
   return m_ca;
}

/*----------------------------------------------------------------------------*/

Atom* Residue::ca() {
   return m_ca;
}

/*----------------------------------------------------------------------------*/

const Atom* Residue::cb() const {
   return m_cb;
}

/*----------------------------------------------------------------------------*/

Atom* Residue::cb() {
   return m_cb;
}

/*----------------------------------------------------------------------------*/

const Atom* Residue::c() const {
   return m_c;
}

/*----------------------------------------------------------------------------*/

Atom* Residue::c() {
   return m_c;
}

/*----------------------------------------------------------------------------*/

const Atom* Residue::n() const {
   return m_n;
}

/*----------------------------------------------------------------------------*/

Atom* Residue::n() {
   return m_n;
}

/*----------------------------------------------------------------------------*/

const Atom* Residue::o() const {
   return m_o;
}

/*----------------------------------------------------------------------------*/

Atom* Residue::o() {
   return m_o;
}

/*----------------------------------------------------------------------------*/

std::string Residue::toString() const {
   return common::s_printf("%c %s%c %s", atoms()[0]->chainId(),
         atoms()[0]->resSeq(), atoms()[0]->chainId(), atoms()[0]->resName());
}

/*----------------------------------------------------------------------------*/

static inline void addDistChecked(std::vector<Residue::Bond*> *bonds, Atom *a1, Atom *a2) {
   static const float maxbonddist_sq = 2.56f; /* 1.6^2 */

   if (a1->squareDist(*a2) < maxbonddist_sq) {
      bonds->push_back(new Residue::Bond(a1, a2));
   }
}

/*----------------------------------------------------------------------------*/

void Residue::calcBonds() {
   for (unsigned int l = 0; l < atoms().size(); l++) {
      for (unsigned int k = l + 1; k < atoms().size(); k++) {
         addDistChecked(&m_bonds, atoms()[k], atoms()[l]);
      }
   }
}

/*----------------------------------------------------------------------------*/

const char* Residue::getName() const {
   return atoms()[0]->resName();
}

/*----------------------------------------------------------------------------*/

Vector Residue::getCbPos() const {
   /* do we have a Cb? */
   if (this->m_cb != NULL) {
      return Vector(*this->m_cb);
   } else if ((this->m_ca != NULL) && (this->m_c != NULL) && (this->m_n != NULL)) {
      Vector ca, c, n;
      Vector cn05, cn05ca;
      Vector cb;

      /* calculate Cb position */
      ca = *this->m_ca;
      c = *this->m_c;
      n = *this->m_n;
      /* pos between C and N */
      cn05 = c + (n - c) / 2;
      // Cb point in plane
      cn05ca = ca - cn05;
      cn05ca = ca + cn05ca;
      /* get normal of plane */
      cb = Vector::crossProd((ca - n), (ca - c));
      /* scale to empirical value */
      cb.normalise();
      cb *= 1.215095;
      /* get cb */
      cb = cn05ca + cb;
      return cb;
   } else {
      return Vector(0, 0, 0);
   }
}

/*----------------------------------------------------------------------------*/

char Residue::sseType() const {
   return m_ssetype;
}

/*----------------------------------------------------------------------------*/

char & Residue::sseType() {
   return m_ssetype;
}

/*----------------------------------------------------------------------------*/

std::vector<Residue::Bond*> Residue::bonds() const {
    return m_bonds;
}

/*----------------------------------------------------------------------------*/

std::vector<Residue::Bond*> & Residue::bonds() {
    return m_bonds;
}

/*----------------------------------------------------------------------------*/

void Residue::destroy(std::vector<Residue*> *resis) {
   for (unsigned int i = 0; i < resis->size(); i++) {
      delete (*resis)[i];
   }
   resis->clear();
}

/*----------------------------------------------------------------------------*/

Vector Residue::caCentre(const Residues &resis) {
   Vector centre(0, 0, 0);
   unsigned int num;

   num = 0;
   for (unsigned int i = 0; i < resis.size(); i++) {
      if (resis[i]->ca() != NULL) {
         centre += *resis[i]->ca();
         num++;
      }
   }
   if (num != 0) {
      centre /= num;
   }
   return centre;
}

/*----------------------------------------------------------------------------*/

Vector Residue::centre(const Residues &resis) {
   Vector centre(0, 0, 0);
   unsigned int num;

   num = 0;
   for (unsigned int i = 0; i < resis.size(); i++) {
      for (unsigned int j = 0; j < resis[i]->atoms().size(); j++) {
         centre += *resis[i]->atoms()[j];
         num++;
      }
   }
   if (num != 0) {
      centre /= num;
      return centre;
   }
   return Vector(0, 0, 0);
}


/*----------------------------------------------------------------------------*/

bool Residue::hasContact(const Residue *res1, const Residue *res2,
      float maxdist) {
   /* we dont want to solve the sqrt */
   const float dist = maxdist * maxdist;
   for (unsigned int i = 0; i < res1->atoms().size(); i++) {
      for (unsigned int j = 0; j < res2->atoms().size(); j++) {
         if (res1->atoms()[i]->squareDist(*res2->atoms()[j]) < dist) {
            return true;
         }
      }
   }
   return false;
}

/*----------------------------------------------------------------------------*/

bool Residue::hasContact_fast(const Residue *res1, const Residue *res2,
      float maxdist) {

#define MAX_CONTACT_CA_DIST_SQ 289

   /* require ca atoms to be close enough, before we look at the others */
   if (res1->hasCa() && res2->hasCa()) {
      if (res1->ca()->squareDist(*res2->ca()) > MAX_CONTACT_CA_DIST_SQ) {
         return false;
      }
   }

   /* we dont want to solve the sqrt */
   const float dist = maxdist * maxdist;
   for (unsigned int i = 0; i < res1->atoms().size(); i++) {
      for (unsigned int j = 0; j < res2->atoms().size(); j++) {
         if (res1->atoms()[i]->squareDist(*res2->atoms()[j]) < dist) {
            return true;
         }
      }
   }
   return false;
}

/*----------------------------------------------------------------------------*/

std::string Residue::getResSequence(const Residues &resis) {
   std::string sequence(resis.size(), '\0');

   for (unsigned int i = 0; i < resis.size(); i++) {
      if (resis[i]->isAminoacid() == true) {
         sequence[i] = AaMap::getResn1(resis[i]->atoms()[0]->resName());
      }
#define ALLOW_GAPS
      else {
#ifdef ALLOW_GAPS
         /* insert gaps for non amino acid residues */
         sequence[i] = '_';
#else
         sequence[i] = AaMap::getResn1(resis[i]->atoms()[0]->resName());
#endif
      }
   }
   return sequence;
}

/*----------------------------------------------------------------------------*/

float Residue::getRmsd(std::vector<Residue*> *m1, std::vector<Residue*> *m2) {
   unsigned int minres;
   float rmsd;

   rmsd = 0;
   minres = std::min((*m1).size(), (*m2).size());
   for (unsigned int i = 0; i < minres; i++) {
      assert(((*m1)[i]->m_ca != NULL) && ((*m2)[i]->m_ca != NULL));
      rmsd += (*m1)[i]->m_ca->squareDist(*(*m2)[i]->m_ca);
   }
   return sqrt(rmsd / minres);
}

/*----------------------------------------------------------------------------*/

bool Residue::hasCa() const {
   return (m_ca != NULL);
}

/*----------------------------------------------------------------------------*/

bool Residue::hasCompleteBackbone() const {
   return ((m_ca != NULL) && (m_c != NULL) && (m_o != NULL) && (m_n != NULL));
}

/*----------------------------------------------------------------------------*/

bool Residue::isAminoacid() const {
   return atoms()[0]->isAminoacid();
}

/*----------------------------------------------------------------------------*/

bool Residue::isNucleotide() const {
   return atoms()[0]->isNuclotide();
}

/*----------------------------------------------------------------------------*/

bool Residue::hasStandardCa() const {
   return m_ca != NULL && m_ca->isHetatm() == false;
}

/*----------------------------------------------------------------------------*/

void Residue::getCaResidues(const std::vector<Residue*> &resis, std::vector<
      Residue*> &caresis) {
   caresis.clear();
   for (unsigned int i = 0; i < resis.size(); i++) {
      if (resis[i]->hasCa() == true) {
         caresis.push_back(resis[i]);
      }
   }
}
