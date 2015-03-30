#ifndef RESALIGNER_H_
#define RESALIGNER_H_

#include "../libxtaldata/residue.h"
#include <stdint.h>

/*----------------------------------------------------------------------------*/

typedef char cbool;

/*----------------------------------------------------------------------------*/

class ResAligner {

public:

   struct Alignment {
      std::vector<Residue*> al1;
      std::vector<Residue*> al2;

      std::vector<Residue*> alfull1;
      std::vector<Residue*> alfull2;
      std::vector<cbool>    match;
      uint32_t dist;
   };



   /** @brief calculates an alignment of two sequences using unit costs
    * @param seq1 sequence 1
    * @param seq2 sequence 2
    * @alignment result is stored in alignment[0] and alignment[1]
    * @return distance of the two sequences */
   static unsigned int align(const std::vector<Residue*> &p1,
         const std::vector<Residue*> &p2, std::vector<Residue*> &al1,
         std::vector<Residue*> &al2);


   static Alignment align(const std::vector<Residue*> &p1,
         const std::vector<Residue*> &p2);


protected:
	   enum EdgeType {
	       EDGETYPE_INS,
	       EDGETYPE_DEL,
	       EDGETYPE_REP,
	       EDGETYPE_COUNT
	    };

    const std::vector<Residue*> &m_seq1;
    const std::vector<Residue*> &m_seq2;
    uint32_t m_n1, m_m1; /* length of sequences + 1*/
    Residue* m_gapChar; /* stands for a gap; is not allowed to occur in sequences */

    /* highest value of matrix */
    uint32_t m_imax;
    uint32_t m_jmax;
    int32_t m_smax;

    int32_t * __restrict__ m_dptable; /* costs table */
    EdgeType * __restrict__ m_betable; /* backward edges */

    ResAligner(const std::vector<Residue*> &prot1,
          const std::vector<Residue*> &prot2);
    ~ResAligner();
    unsigned int fillTable();
    void printTables();
    void goBackwards(std::vector<Residue*> &al1, std::vector<Residue*> &al2);
    void goBackwards(Alignment &al);
};

#endif /* RESALIGNER_H_ */
