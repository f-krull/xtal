#include "iresidues.h"
#include <algorithm>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include "matrix.cpp"

/*----------------------------------------------------------------------------*/

IResidues::IResidues() {
}

/*----------------------------------------------------------------------------*/

bool IResidues::build(IAtoms &a) {
   std::string lastres;
   char lastICode;
   char lastchain;
   const char *lastResn = NULL;
   Residue *res;

   /* algorithm needs at least one atom to work */
   if (a.size() == 0) {
      return false;
   }

   /* initialize with first atom */
   std::vector<Atom*>::iterator it_start = a.atoms().begin();
   lastchain = (*it_start)->chainId();
   lastres = (*it_start)->resSeq();
   lastICode =  (*it_start)->iCode();
   lastResn = (*it_start)->resName();

   std::vector<Atom*>::iterator it_curr = a.atoms().begin();
   while (it_curr != a.atoms().end()) {
      /* go to next atom */
      it_curr++;
      /* see if we have a residue completed */

      bool newRes = false;
      newRes = newRes || it_curr == a.atoms().end();
      newRes = newRes || lastchain != (*it_curr)->chainId();
      newRes = newRes || lastres != (*it_curr)->resSeq();
      newRes = newRes || lastICode != (*it_curr)->iCode();
      newRes = newRes || strcmp(lastResn, (*it_curr)->resName()) != 0;

      if (newRes == true) {
         /* pass all belonging atoms */
         res = new Residue(it_start, it_curr);
         /* sort out residues with no residue name */
         if (res->isAminoacid()) {
            resis().push_back(res);
         } else {
            nonResis().push_back(res);
         }
         /* set new start */
         it_start = it_curr;
         if (it_start != a.atoms().end()) {
            lastres = (*it_start)->resSeq();
            lastResn = (*it_start)->resName();
            lastICode = (*it_start)->iCode();
            lastchain = (*it_start)->chainId();
         }
      }
   }
   return true;
}

/*----------------------------------------------------------------------------*/

static inline void deleteResidue(Residue *s) {
   delete s;
}

/*----------------------------------------------------------------------------*/

Vector IResidues::centre() const {
   return Residue::caCentre(resis());
}

/*----------------------------------------------------------------------------*/

void IResidues::destroy() {
   std::for_each(m_resis.begin(), m_resis.end(), deleteResidue);
   std::for_each(m_nonRes.begin(), m_nonRes.end(), deleteResidue);
   m_resis.clear();
   m_nonRes.clear();
   update();
}

/*----------------------------------------------------------------------------*/

IResidues::~IResidues() {
}

/*----------------------------------------------------------------------------*/

const Residues & IResidues::resis() const {
   return m_resis;
}

/*----------------------------------------------------------------------------*/

Residues & IResidues::resis() {
   return m_resis;
}

/*----------------------------------------------------------------------------*/

const Residues & IResidues::nonResis() const {
   return m_nonRes;
}

/*----------------------------------------------------------------------------*/

Residues & IResidues::nonResis() {
   return m_nonRes;
}

/*----------------------------------------------------------------------------*/

size_t IResidues::size() const {
   return m_resis.size();
}

/*----------------------------------------------------------------------------*/

void IResidues::add(Residue *a) {
   m_resis.push_back(a);
   update();
}

/*----------------------------------------------------------------------------*/

void IResidues::update() {
}

/*----------------------------------------------------------------------------*/

void IResidues::calcBonds() {
   const float maxbonddist_sq = 2.56f; /* 1.6^2 */

   for (unsigned int i = 0; i < resis().size(); i++) {
      resis()[i]->calcBonds();
   }
   for (unsigned int i = 0; i + 1 < resis().size(); i++) {
      if ((resis()[i]->c() == NULL) || (resis()[i + 1]->n() == NULL)) {
         continue;
      }
      if (resis()[i]->c()->squareDist(*resis()[i + 1]->n()) < maxbonddist_sq) {
         resis()[i]->bonds().push_back(new Residue::Bond(resis()[i]->c(), resis()[i + 1]->n()));
      }
   }
}

/*----------------------------------------------------------------------------*/

inline float getHbondEnergy(const Vector &c, const Vector &o, const Vector &n,
      const Vector &h) {
   static const float q1 = 0.42;
   static const float q2 = 0.20;
   static const float f = 332;

   return q1 * q2 * (1 / o.dist(n) + 1 / c.dist(h) - 1 / o.dist(h) - 1
         / c.dist(n)) * f;
}

/*----------------------------------------------------------------------------*/

void IResidues::calcSses() {
   typedef char cbool;
   Matrix<bool> hbond(resis().size(), resis().size());
   std::vector<Atom*> hbonds;
   Residue *resi, *resj;
   std::vector<Vector> nhpos;
   std::vector<Residue*> caresis;

   /* select residues with backbone only */
   caresis.reserve(resis().size());
   for (unsigned int i = 0; i < resis().size(); i++) {
      if (resis()[i]->hasCompleteBackbone() == true) {
         caresis.push_back(resis()[i]);
      }
   }
   /* compute N hydrogens */
   Vector c, n, ca, h;
   nhpos.reserve(caresis.size());
   if (caresis.size() > 0) {
      nhpos.push_back(*caresis[0]->n());
   }
   for (unsigned int i = 1; i < caresis.size(); i++) {
      c = *caresis[i - 1]->c();
      n = *caresis[i]->n();
      ca = *caresis[i]->ca();
      h = c + (ca - c) / 2;
      h = n - h;
      h.normalise();
      h = n + h;
      nhpos.push_back(h);
   }
   for (unsigned int i = 0; i < caresis.size(); i++) {
      for (unsigned int j = 0; j < caresis.size(); j++) {
         if (abs((int) i - (int) j) < 3) {
            continue;
         }
         resi = caresis[i];
         resj = caresis[j];
         if (resi->ca()->dist(*resj->ca()) > 10) {
            /* min ~3.5A; max ~8.5A */
            continue;
         }
         if (getHbondEnergy(*resi->c(), *resi->o(), *resj->n(), nhpos[j])
               < -0.5f) {
            hbond[i][j] = true;

         }
      }
   }
   std::vector<cbool> turn4(caresis.size(), false);
   //   vector<cbool> turn3((*resis).size(), false);
   //   vector<cbool> turn5((*resis).size(), false);
   //   for (unsigned int i = 0; i + 3 < (*resis).size(); i++) {
   //      if (hbond[i][i + 3] == true) {
   //         turn3[i] = true;
   //      }
   //   }
   for (unsigned int i = 0; i + 4 < caresis.size(); i++) {
      if (hbond[i][i + 4] == true) {
         turn4[i] = true;
      }
   }
   //   for (unsigned int i = 0; i + 5 < (*resis).size(); i++) {
   //      if (hbond[i][i + 5] == true) {
   //         turn5[i] = true;
   //      }
   //   }
   for (unsigned int i = 0; i < caresis.size(); i++) {
      caresis[i]->sseType() = Residue::COIL;
   }
   //   for (unsigned int i = 0; i + 3 < (*resis).size(); i++) {
   //      if (turn3[i] && turn3[i + 3]) {
   //         (*resis)[i + 0]->ssetype = Residue::HELIX;
   //         (*resis)[i + 1]->ssetype = Residue::HELIX;
   //         (*resis)[i + 2]->ssetype = Residue::HELIX;
   //         (*resis)[i + 3]->ssetype = Residue::HELIX;
   //      }
   //   }
   for (unsigned int i = 0; i + 4 < caresis.size(); i++) {
      if (turn4[i] && turn4[i + 4]) {
         caresis[i + 0]->sseType() = Residue::HELIX;
         caresis[i + 1]->sseType() = Residue::HELIX;
         caresis[i + 2]->sseType() = Residue::HELIX;
         caresis[i + 3]->sseType() = Residue::HELIX;
         caresis[i + 4]->sseType() = Residue::HELIX;
      }
   }
   //   for (unsigned int i = 0; i + 5 < (*resis).size(); i++) {
   //      if (turn5[i] && turn5[i + 5]) {
   //         (*resis)[i + 0]->ssetype = Residue::HELIX;
   //         (*resis)[i + 1]->ssetype = Residue::HELIX;
   //         (*resis)[i + 2]->ssetype = Residue::HELIX;
   //         (*resis)[i + 3]->ssetype = Residue::HELIX;
   //         (*resis)[i + 4]->ssetype = Residue::HELIX;
   //         (*resis)[i + 5]->ssetype = Residue::HELIX;
   //      }
   //   }
   //   Matrix<bool> pbridge(resis->size(), resis->size(), false);
   //   Matrix<bool> abridge(resis->size(), resis->size(), false);
   for (unsigned int i = 1; i + 1 < caresis.size(); i++) {
      for (unsigned int j = 1; j + 1 < caresis.size(); j++) {
         /* parallel bridge */
         if (hbond[i - 1][j] && hbond[j][i + 1]) {
            //pbridge[i][j] = true;
            caresis[i]->sseType() = Residue::SHEET;
            caresis[j]->sseType() = Residue::SHEET;
         } else if (hbond[j - 1][i] && hbond[i][j + 1]) {
            //pbridge[i][j] = true;
            caresis[i]->sseType() = Residue::SHEET;
            caresis[j]->sseType() = Residue::SHEET;
         }
         /* antiparallel bridge */
         if (hbond[i][j] && hbond[j][i]) {
            //abridge[i][j] = true;
            caresis[i]->sseType() = Residue::SHEET;
            caresis[j]->sseType() = Residue::SHEET;
         } else if (hbond[i - 1][j + 1] && hbond[j - 1][i + 1]) {
            //abridge[i][j] = true;
            caresis[i]->sseType() = Residue::SHEET;
            caresis[j]->sseType() = Residue::SHEET;
         }
      }
   }
   for (unsigned int i = 1; i + 1 < caresis.size(); i++) {
      if (caresis[i]->sseType() == Residue::SHEET) {
         if ((caresis[i - 1]->sseType() != Residue::SHEET)
               && (caresis[i + 1]->sseType() != Residue::SHEET)) {
            caresis[i]->sseType() = Residue::COIL;
         }
      }
   }
}


