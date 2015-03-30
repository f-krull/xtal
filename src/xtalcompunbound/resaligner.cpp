#include "resaligner.h"
#include <algorithm>
#include <string.h>
#include <assert.h>

using namespace std;

/*----------------------------------------------------------------------------*/

ResAligner::ResAligner(const std::vector<Residue*> &p1,
  const std::vector<Residue*> &p2) : m_seq1(p1), m_seq2(p2) {
   m_n1 = p1.size() + 1;
   m_m1 = p2.size() + 1;
   m_dptable = new int32_t[m_n1 * m_m1];
   memset(m_dptable, 0, sizeof(int32_t) * m_n1 * m_m1);
   m_betable = new EdgeType[m_n1 * m_m1];
   memset(m_betable, 0, sizeof(EdgeType) * m_n1 * m_m1);
   this->m_gapChar = NULL;
   this->m_imax = 0;
   this->m_jmax = 0;
   this->m_smax = 0;
}

/*----------------------------------------------------------------------------*/

ResAligner::~ResAligner() {
   delete [] m_dptable;
   delete [] m_betable;
}

/*----------------------------------------------------------------------------*/

static int32_t getCost(const Residue* res1, const Residue *res2) {
   if (res2 == NULL || res1 == NULL) {
      return -1;
   }
   /* similar names? */
   if (res1->atoms()[0]->isResn(res2->atoms()[0]->resName())) {
	   return 2;
   }
   return -1;
}

/*----------------------------------------------------------------------------*/

uint32_t ResAligner::fillTable() {
	int32_t tmp_costs;
	   int32_t* __restrict__ dpcurr;
	   EdgeType* __restrict__ becurr;

	   assert (m_gapChar != m_seq1[0]);
	   assert (m_gapChar != m_seq2[0]);
	   m_dptable[0] = 0;
	   m_betable[0] = EDGETYPE_REP;
	   for (uint32_t i = 1; i < m_n1; i++) {
	      /* gap in seq2 */
	      m_dptable[i * m_m1] = m_dptable[(i - 1) * m_m1] + getCost(m_seq1[i - 1], m_gapChar);
	      m_betable[i * m_m1] = EDGETYPE_DEL;
	      /* start of new alignment? */
	      m_dptable[i * m_m1] = std::max(m_dptable[i * m_m1], 0);
	      /* end of local alignment */
	      if (m_dptable[i * m_m1] > m_smax) {
	         m_smax = m_dptable[i * m_m1];
	         m_imax = i;
	         m_jmax = 0;
	      }
	   }
	   for (uint32_t j = 1; j < m_m1; j++) {
	      /* gap in seq1 */
	      m_dptable[0 + j] = m_dptable[0 + (j - 1)] + getCost(m_gapChar, m_seq2[j - 1]);
	      m_betable[0 + j] = EDGETYPE_INS;
	      /* start of new alignment? */
	      m_dptable[0 + j] = std::max(m_dptable[0 + j], 0);
	      if (m_dptable[0 + j] > m_smax) {
	         m_smax = m_dptable[0 + j];
	         m_imax = 0;
	         m_jmax = j;
	      }

	      for (uint32_t i = 1; i < m_n1; i++) {
	         /* short cut to position of matrix we want to fill */
	         dpcurr = &m_dptable[i * m_m1 + j];
	         becurr = &m_betable[i * m_m1 + j];

	         /* gap in seq1 */
	         tmp_costs = m_dptable[i * m_m1 + (j - 1)]
	               + getCost(m_gapChar, m_seq2[j - 1]);
	         if (tmp_costs >= dpcurr[0]) {
	            dpcurr[0] = tmp_costs;
	            becurr[0] = EDGETYPE_INS;
	         }
	         /* gap in seq2 */
	         tmp_costs = m_dptable[(i - 1) * m_m1 + j]
	               + getCost(m_seq1[i - 1], m_gapChar);
	         if (tmp_costs >= dpcurr[0]) {
	            dpcurr[0] = tmp_costs;
	            becurr[0] = EDGETYPE_DEL;
	         }
	         /* replacement */
	         tmp_costs = m_dptable[(i - 1) * m_m1 + (j - 1)]
	               + getCost(m_seq1[i - 1], m_seq2[j - 1]);
	         if (tmp_costs >= dpcurr[0]) {
	            dpcurr[0] = tmp_costs;
	            becurr[0] = EDGETYPE_REP;
	         }
	         /* end of local alignment? */
	         if (dpcurr[0] > m_smax) {
	            m_smax = dpcurr[0];
	            m_imax = i;
	            m_jmax = j;
	         }
	      }
	   }
	   return m_dptable[m_imax * m_m1 + m_jmax];
}

/*----------------------------------------------------------------------------*/

void ResAligner::goBackwards(std::vector<Residue*> &al1,
		std::vector<Residue*> &al2) {
	uint32_t i, j;

	al1.clear();
	al2.clear();
	i = m_imax;
	j = m_jmax;
	/* which path was a good alignment? (greedy - multiple choices) */
	while ((i > 0) || (j > 0)) {
		if (m_dptable[i * m_m1 + j] == 0) {
			/* stop if begin of local alignment is reached */
			break;
		} else if (m_betable[i * m_m1 + j] == EDGETYPE_INS) {
				j--;
			} else if (m_betable[i * m_m1 + j] == EDGETYPE_DEL) {
				i--;
			} else { /* betable[i][j] == REP */
				if (m_seq1[i - 1]->atoms()[0]->isResn(
						m_seq2[j - 1]->atoms()[0]->resName()) == true) {
					al1.push_back(m_seq1[i - 1]);
					al2.push_back(m_seq2[j - 1]);
				}
				i--;
				j--;
			}
   }
   /* reverse result */
   std::reverse(al1.begin(), al1.end());
   std::reverse(al2.begin(), al2.end());
}

/*----------------------------------------------------------------------------*/

void ResAligner::goBackwards(ResAligner::Alignment &al) {
   uint32_t i, j;

   al.al1.clear();
   al.al2.clear();
   al.alfull1.clear();
   al.alfull2.clear();
   al.match.clear();

   i = (m_n1 > 0 ? m_n1 -1 : 0);
   j = (m_m1 > 0 ? m_m1 -1 : 0);


   while ((i > m_imax) || (j > m_jmax)) {
      /* i-- */
      if (i - m_imax < j - m_jmax) {
         al.alfull1.push_back(NULL);
         al.alfull2.push_back(m_seq2[j - 1]);
         al.match.push_back(false);
         j--;
      } else if (i - m_imax > j - m_jmax) {
         al.alfull1.push_back(m_seq1[i - 1]);
         al.alfull2.push_back(NULL);
         al.match.push_back(false);
         i--;
      } else {
         al.alfull1.push_back(m_seq1[i - 1]);
         al.alfull2.push_back(m_seq2[j - 1]);
         al.match.push_back(false);
         i--;
         j--;
      }
   }



   /* which path was a good alignment? (greedy - multiple choices) */
   while ((i > 0) || (j > 0)) {
      if (m_dptable[i * m_m1 + j] == 0) {
         /* stop if begin of local alignment is reached */
         break;
      } else if (m_betable[i * m_m1 + j] == EDGETYPE_INS) {
            al.alfull1.push_back(NULL);
            al.alfull2.push_back(m_seq2[j - 1]);
            al.match.push_back(false);
            j--;
         } else if (m_betable[i * m_m1 + j] == EDGETYPE_DEL) {
            al.alfull1.push_back(m_seq1[i - 1]);
            al.alfull2.push_back(NULL);
            al.match.push_back(false);
            i--;
         } else { /* betable[i][j] == REP */
            if (m_seq1[i - 1]->atoms()[0]->isResn(
                  m_seq2[j - 1]->atoms()[0]->resName()) == true) {
               al.al1.push_back(m_seq1[i - 1]);
               al.al2.push_back(m_seq2[j - 1]);
               al.match.push_back(true);
            } else {
               al.match.push_back(false);
            }
            al.alfull1.push_back(m_seq1[i - 1]);
            al.alfull2.push_back(m_seq2[j - 1]);
            i--;
            j--;
         }
   }

   while ((i > 0) || (j > 0)) {
      if (i == 0 && j > 0) {
         al.alfull1.push_back(NULL);
         al.alfull2.push_back(m_seq2[j - 1]);
         al.match.push_back(false);
         j--;
      } else if (i > 0 && j == 0) {
         al.alfull1.push_back(m_seq1[i - 1]);
         al.alfull2.push_back(NULL);
         al.match.push_back(false);
         i--;
      } else {
         al.alfull1.push_back(m_seq1[i - 1]);
         al.alfull2.push_back(m_seq2[j - 1]);
         al.match.push_back(false);
         i--;
         j--;
      }
   }


   /* reverse result */
   std::reverse(al.al1.begin(), al.al1.end());
   std::reverse(al.al2.begin(), al.al2.end());
   std::reverse(al.alfull1.begin(), al.alfull1.end());
   std::reverse(al.alfull2.begin(), al.alfull2.end());
   std::reverse(al.match.begin(), al.match.end());
}

/*----------------------------------------------------------------------------*/

uint32_t ResAligner::align(const std::vector<Residue*> &p1,
      const std::vector<Residue*> &p2, std::vector<Residue*> &al1,
      std::vector<Residue*> &al2) {
   uint32_t dist;

   ResAligner ra(p1, p2);
   dist = ra.fillTable();
   ra.goBackwards(al1, al2);
   return dist;
}

/*----------------------------------------------------------------------------*/

ResAligner::Alignment ResAligner::align(const std::vector<Residue*> &p1,
      const std::vector<Residue*> &p2) {
   Alignment al;

   ResAligner ra(p1, p2);
   al.dist = ra.fillTable();
   ra.goBackwards(al);
   return al;
}
