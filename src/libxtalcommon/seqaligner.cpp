#include "seqaligner.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string.h>
#include <assert.h>
#include <new>

/*----------------------------------------------------------------------------*/

SeqAligner::SeqAligner(const char *seq1, const char *seq2,
                       char gapChar) {
   m_n1 = strlen(seq1) + 1;
   m_m1 = strlen(seq2) + 1;
   try {
      m_dptable = new int32_t[m_n1 * m_m1];
      memset(m_dptable, 0, sizeof(int32_t) * m_n1 * m_m1);
      m_betable = new EdgeType[m_n1 * m_m1];
      memset(m_betable, 0, sizeof(EdgeType) * m_n1 * m_m1);
   } catch (std::bad_alloc& ba)
   {
     std::cout << "bad_alloc caught: " << ba.what() << std::endl;
     std::cout << "seq1 size " << m_n1 << " : " << seq1 << std::endl;
     std::cout << "seq2 size " << m_m1 << " : " << seq2 << std::endl;
   }

   this->m_seq1 = seq1;
   this->m_seq2 = seq2;
   this->m_gapChar = gapChar;
   this->m_imax = 0;
   this->m_jmax = 0;
   this->m_smax = 0;
}

/*----------------------------------------------------------------------------*/

SeqAligner::~SeqAligner() {
   delete [] m_dptable;
   delete [] m_betable;
}

/*----------------------------------------------------------------------------*/

void SeqAligner::printTables() {
	std::cout << "   " << "  " << m_gapChar;
	for (uint32_t i = 0; i+1 < m_n1; i++) {
		std::cout << std::setw(3) << m_seq1[i];
	}
	std::cout << std::endl;
	std::cout << "  " << m_gapChar;
	for (uint32_t i = 0; i < m_n1; i++) {
		std::cout << std::setw(3) << m_dptable[i * m_m1 + 0];
	}
	std::cout << std::endl;
	for (uint32_t j = 1; j < m_m1; j++) {
		std::cout << std::setw(3) << m_seq2[j - 1];
		for (uint32_t i = 0; i < m_n1; i++) {
			std::cout << std::setw(3) << m_dptable[i * m_m1 + j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << "   " << "  " << m_gapChar;
	for (uint32_t i = 0; (i+1) < m_n1; i++) {
		std::cout << std::setw(3) << m_seq1[i];
	}

	/* map enum to chars */
	char egdeName[EDGETYPE_COUNT];
	egdeName[EDGETYPE_INS] = 'i';
	egdeName[EDGETYPE_DEL] = 'd';
	egdeName[EDGETYPE_REP] = 'r';

	std::cout << std::endl;
	std::cout << "  " << m_gapChar;
	for (uint32_t i = 0; i < m_n1; i++) {
		std::cout << std::setw(3) << egdeName[m_betable[i * m_m1 + 0]];
	}
	std::cout << std::endl;
	for (uint32_t j = 1; j < m_m1; j++) {
		std::cout << std::setw(3) << m_seq2[j - 1];
		for (uint32_t i = 0; i < m_n1; i++) {
			std::cout << std::setw(3) << egdeName[m_betable[i * m_m1 + j]];
		}
		std::cout << std::endl;
	}
}

/*----------------------------------------------------------------------------*/

void SeqAligner::goBackwards(std::pair<std::string, std::string> &alignment) {
	char *s1, *s2;
	uint32_t i, j, k;

	i = m_imax;
	j = m_jmax;
	k = 0;
	s1 = new char[m_n1 + m_m1];
	s2 = new char[m_n1 + m_m1];
	/* which path was a good alignment? (greedy - multiple choices) */


	while ((i > 0) || (j > 0)) {
	   if (m_dptable[i * m_m1 + j] == 0) {
	      /* stop if begin of local alignment is reached */
	      break;
	   } else if (m_betable[i * m_m1 + j] == EDGETYPE_INS) {
			s2[k] = (m_seq2[j - 1]);
			s1[k] = m_gapChar;
			j--;
		} else if (m_betable[i * m_m1 + j] == EDGETYPE_DEL) {
			s2[k] = m_gapChar;
			s1[k] = m_seq1[i - 1];
			i--;
		} else { /* betable[i][j] == REP */
			s1[k] = m_seq1[i - 1];
			s2[k] = m_seq2[j - 1];
			i--;
			j--;
		}
		k++;
	}

	/* terminate string */
	s1[k] = '\0';
	s2[k] = '\0';
	alignment.second = s1;
	alignment.first = s2;
	reverse(alignment.second.begin(), alignment.second.end());
	reverse(alignment.first.begin(), alignment.first.end());
	delete[] s1;
	delete[] s2;
}



/*----------------------------------------------------------------------------*/
//
//template <typename _Costs>
//float SeqAligner<_Costs>::getSeqId(const std::string &seq1, const std::string &seq2, char gapChar) {
//   int32_t max, min;
//	max = std::max(seq1.size(), seq2.size());
//	min = std::min(seq1.size(), seq2.size());
//	if (min == 0) {
//		return 0;
//	}
//	return ((float) max - getScore(seq1.data(), seq2.data(), gapChar)) / min;
//}

/*----------------------------------------------------------------------------*/

