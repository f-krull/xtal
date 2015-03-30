#ifndef RESIDUE_H_
#define RESIDUE_H_

#include "iatoms.h"

/*----------------------------------------------------------------------------*/

class Residue;

/*----------------------------------------------------------------------------*/

typedef std::vector<Residue*> Residues;

/*----------------------------------------------------------------------------*/

class Residue: public IAtoms {
public:

   class Bond {
   private:
   protected:
   public:
      Atom* atm1;
      Atom* atm2;

      Bond(Atom *a, Atom *b);
   };


   static char HELIX;
   static char SHEET;
   static char COIL;
   static char UNSET;

   Residue();
   Residue(const Residue &r2);
   Residue(Atoms::iterator &b, Atoms::iterator &end);
   virtual ~Residue();



   const char* getName() const;
   std::string toString() const;

   Atom* ca();
   Atom* cb();
   Atom* c();
   Atom* o();
   Atom* n();
   const Atom* ca() const;
   const Atom* cb() const;
   const Atom* c() const;
   const Atom* o() const;
   const Atom* n() const;
   Vector getCbPos() const;



   char & sseType();
   char sseType() const;

   std::vector<Bond*> bonds() const;
   std::vector<Bond*> & bonds();

   bool hasCa() const;
   bool hasCompleteBackbone() const;
   bool isAminoacid() const;
   bool isNucleotide() const;
   bool hasStandardCa() const; /* Ca is not HETATM */

   Residue & operator=(const Residue &p2);

   void calcBonds();


   static void destroy(Residues *resis);
   static void destroyDeep(Residues *resis);
   static Vector centre(const Residues &resis);
   static Vector caCentre(const Residues &resis);
   static bool hasContact(const Residue *res1, const Residue *res2,
         float maxdist = 3.5);
   static bool hasContact_fast(const Residue *res1, const Residue *res2,
         float maxdist = 3.5);
   static std::string getResSequence(const Residues &resis);
   static float getRmsd(Residues *m1, Residues *m2);
   static void getCaResidues(const Residues &resis, Residues &caresis);
protected:
   void init();

   Atom *m_ca;
   Atom *m_c;
   Atom *m_o;
   Atom *m_n;
   Atom *m_cb;

   std::vector<Bond*> m_bonds;
   char m_ssetype; /* 'H','E' or 'C' */

};

#endif /* RESIDUE_H_ */
