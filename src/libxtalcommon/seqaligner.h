#ifndef SEQALIGNER_H_
#define SEQALIGNER_H_

#include <string>
#include <stdint.h>
#include <assert.h>

//#define XTALDEBUG

/*----------------------------------------------------------------------------*/

class UnitScore {
public:
   inline int32_t operator()(const char a, const char b) const {
      if (a == b) {
         return 2;
      }
      return -1;
   }
   inline int32_t getBestScore() const {
	   return 2;
   }
};

#include <limits.h>
class UnitScoreCache {
public:
   UnitScoreCache() : m_c() {
      for (uint8_t i = 0; i < UCHAR_MAX; i++) {
         for (uint8_t j = i; j < UCHAR_MAX; j++) {
            m_t[i * UCHAR_MAX + j] = m_t[j * UCHAR_MAX + i] = m_c(i ,j);
         }
      }
   }

   inline int32_t operator()(const uint8_t a, const uint8_t b) const {
      return m_t[a * UCHAR_MAX + b];
   }

   inline int32_t getBestScore() const {
      return m_c.getBestScore();
   }

private:
   int32_t m_t[UCHAR_MAX * UCHAR_MAX];
   const UnitScore m_c;
};

/*----------------------------------------------------------------------------*/

class SeqAligner {
public:

   enum SpeedUp {
      SPEEDUP_NONE,
      SPEEDUP_SSE_16,
      SPEEDUP_SSE_16_AUTO,
      SPEEDUP_SSE_32,
      SPEEDUP_SSE_32_AUTO,
      SPEEDUP_COUNT
   };

   template<typename _Costs>
   static int32_t align(const char *seq1, const char *seq2,
         std::pair<std::string, std::string> *alignment, const _Costs &c
         , SpeedUp s = SPEEDUP_SSE_32, char gapChar = '-');

   template<typename _Costs>
   static float getSeqSim(const std::string &seq1, const std::string &seq2,
         const _Costs &c, SpeedUp s = SPEEDUP_NONE, char gapChar = '-');

protected:
   enum EdgeType {
       EDGETYPE_INS,
       EDGETYPE_DEL,
       EDGETYPE_REP,
       EDGETYPE_COUNT
    };

    SeqAligner(const char *seq1, const char *seq2,
                char gapChar);
    virtual ~SeqAligner();


    template <typename _Costs>
    int32_t fillTable(_Costs costs);
    template <typename _Costs>
    int32_t fillTable_sse_32(_Costs costs);
    template <typename _Costs>
    int32_t fillTable_sse_16(_Costs costs);
    virtual void printTables();
    virtual void goBackwards(std::pair<std::string, std::string> &alignment);

    int32_t * __restrict__ m_dptable; /* costs table */
    EdgeType * __restrict__ m_betable; /* backward edges */

    const char* m_seq1;
    const char* m_seq2;
    uint32_t m_n1, m_m1; /* length of sequences + 1*/
    char m_gapChar; /* stands for a gap; is not allowed to occur in sequences */

    /* highest value of matrix */
    uint32_t m_imax;
    uint32_t m_jmax;
    int32_t m_smax;
};

/*----------------------------------------------------------------------------*/
#include <string.h>
template <typename _Costs>
int32_t SeqAligner::align(const char *seq1, const char *seq2,
         std::pair<std::string, std::string> *alignment, const _Costs &c, SpeedUp s, char gapChar) {
   int32_t score;

   SeqAligner sa(seq1, seq2, gapChar);
#if 0
   if (s == SPEEDUP_SSE_32_AUTO && strlen(seq1) * strlen(seq2) > (300 * 300)) {
      score = sa.fillTable_sse_32(c);
   } else if (s == SPEEDUP_SSE_32) {
      score = sa.fillTable_sse_32(c);
   } else if (s == SPEEDUP_SSE_16_AUTO && strlen(seq1) * strlen(seq2) > (300 * 300)) {
      score = sa.fillTable_sse_16(c);
   } else if (s == SPEEDUP_SSE_16) {
      score = sa.fillTable_sse_16(c);
   } else {
      score = sa.fillTable(c);
   }
#else
   score = sa.fillTable(c);
#endif



#ifdef XTALDEBUG
   sa.printTables();
#endif
   if (alignment != NULL) {
      sa.goBackwards(*alignment);
   }
   return score;
}

/*----------------------------------------------------------------------------*/


template <typename _Costs>
float SeqAligner::getSeqSim(const std::string &seq1, const std::string &seq2, const _Costs &c, SpeedUp s, char gapChar) {
   float score = 0;

   score = SeqAligner::align(seq1.c_str(), seq2.c_str(), NULL, c, s, '-');
   /* normalize by max score */
   score /= c.getBestScore();
   /* normalize by smaller sequence length */
   score /= std::min(seq1.size(), seq2.size());
   return score;
}

/*----------------------------------------------------------------------------*/

template <typename _Costs>
int32_t SeqAligner::fillTable(_Costs costs) {
   int32_t tmp_costs;
   int32_t* __restrict__ dpcurr;
   EdgeType* __restrict__ becurr;

   assert (m_gapChar != m_seq1[0]);
   assert (m_gapChar != m_seq2[0]);
   m_dptable[0] = 0;
   m_betable[0] = EDGETYPE_REP;
   for (uint32_t i = 1; i < m_n1; i++) {
      /* gap in seq2 */
      m_dptable[i * m_m1] = m_dptable[(i - 1) * m_m1] + costs(m_seq1[i - 1], m_gapChar);
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
      m_dptable[0 + j] = m_dptable[0 + (j - 1)] + costs(m_gapChar, m_seq2[j - 1]);
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
               + costs(m_gapChar, m_seq2[j - 1]);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_INS;
         }
         /* gap in seq2 */
         tmp_costs = m_dptable[(i - 1) * m_m1 + j]
               + costs(m_seq1[i - 1], m_gapChar);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_DEL;
         }
         /* replacement */
         tmp_costs = m_dptable[(i - 1) * m_m1 + (j - 1)]
               + costs(m_seq1[i - 1], m_seq2[j - 1]);
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

#ifdef __SSE4_1__
#include <xmmintrin.h>
#include <emmintrin.h>
#include <tmmintrin.h>
#include <pmmintrin.h>

template<typename _Costs>
__attribute__((optimize("unroll-loops")))
int32_t SeqAligner::fillTable_sse_16(_Costs costs ) {
   int32_t tmp_costs;
   int32_t* __restrict__ dpcurr;
   EdgeType* __restrict__ becurr;

   assert(m_gapChar != m_seq1[0]);
   assert(m_gapChar != m_seq2[0]);
   m_dptable[0] = 0;
   m_betable[0] = EDGETYPE_REP;

   const uint32_t V_SIZE = 8u;

   /* init first row ------------------------- */
   for (uint32_t j = 1; j < m_m1; j++) {
      /* gap in seq2 */
      m_dptable[0 + j] = m_dptable[0 + (j - 1)]
            + costs(m_seq1[j - 1], m_gapChar);
      m_betable[0 + j] = EDGETYPE_DEL;
      /* start of new alignment? */
      m_dptable[0 + j] = std::max(m_dptable[0 + j], 0);
      /* end of local alignment */
      if (m_dptable[0 + j] > m_smax) {
         m_smax = m_dptable[0 + j];
         m_imax = 0;
         m_jmax = j;
      }
   }

   /* init first V_SIZE columns-------------------------- */
   for (uint32_t i = 1; i < m_n1; i++) {
      /* gap in seq1 */
      m_dptable[i * m_m1] = m_dptable[(i - 1) * m_m1]
            + costs(m_gapChar, m_seq2[i - 1]);
      m_betable[i * m_m1] = EDGETYPE_INS;
      /* start of new alignment? */
      m_dptable[i * m_m1] = std::max(m_dptable[i * m_m1], 0);
      if (m_dptable[i * m_m1] > m_smax) {
         m_smax = m_dptable[i * m_m1];
         m_imax = 0;
         m_jmax = i;
      }
      for (uint32_t j = 1; j < std::min(m_m1, V_SIZE) + 2; j++) {
         /* short cut to position of matrix we want to fill */
         dpcurr = &m_dptable[i * m_m1 + j];
         becurr = &m_betable[i * m_m1 + j];

         /* gap in seq1 */
         tmp_costs = m_dptable[i * m_m1 + (j - 1)]
               + costs(m_gapChar, m_seq2[j - 1]);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_INS;
         }
         /* gap in seq2 */
         tmp_costs = m_dptable[(i - 1) * m_m1 + j]
               + costs(m_seq1[i - 1], m_gapChar);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_DEL;
         }
         /* replacement */
         tmp_costs = m_dptable[(i - 1) * m_m1 + (j - 1)]
               + costs(m_seq1[i - 1], m_seq2[j - 1]);
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

   /* use diagonal */
#if 1

   union _TD {
      __m128i s;
      int16_t c[V_SIZE];

   } pp, p, t, cd;

   union _TE {
      __m128i s;
      int16_t e[V_SIZE];

   };
   const _TE cbi = {_mm_set1_epi16(EDGETYPE_INS)};
   const _TE cbd = {_mm_set1_epi16(EDGETYPE_DEL)};
   const _TE cbr = {_mm_set1_epi16(EDGETYPE_REP)};

   _TE cb;

   __m128i _b;
   const __m128i  _s0 = _mm_set1_epi16(0);
#define SSECMP16 _mm_cmpgt_epi16


   //char t[V_SIZE];
   uint32_t endi = m_n1 - ((m_n1 - 1) % V_SIZE);
   //endi = std::max(endi, V_SIZE) - V_SIZE;
   uint32_t endj = m_m1 - (V_SIZE);
   for (uint32_t i = 1; i < endi; i += V_SIZE) {
      for (uint32_t j = 2; j < endj; j += 1) {
         if (j == 2) {
            for (uint32_t k = 0; k < V_SIZE; k++) {
               pp.c[k] = m_dptable[(i + (V_SIZE - 1 - k)) * m_m1 + (j + k - 2)];
               p.c[k] = m_dptable[(i + (V_SIZE - 1 - k)) * m_m1 + (j + k - 1)];
            }
         } else {
            {
               pp = p;
            }
            for (uint32_t k = 0; k < V_SIZE; k++) {
               p.c[k] = m_dptable[(i + (V_SIZE - 1 - k)) * m_m1 + (j + k - 1)];
            }
         }
         cb.s = cbi.s;
         cd.s = _s0;

         /* INS */
         for (uint32_t k = 0; k < (V_SIZE - 0); k++) {
            t.c[k] = p.c[k] + costs(m_gapChar, m_seq2[j + k - 1]);
         }
         _b = SSECMP16(t.s, cd.s);
         cd.s = _mm_or_si128(_mm_and_si128(_b, t.s), _mm_andnot_si128(_b, cd.s));
         cb.s = _mm_or_si128(_mm_and_si128(_b, cbi.s), _mm_andnot_si128(_b, cb.s));


         /* DEL */
         for (uint32_t k = 0; k < (V_SIZE - 1); k++) {
            t.c[k] = p.c[k + 1] + costs(m_seq1[i + (V_SIZE - 1 - k) - 1], m_gapChar);
         }
         t.c[V_SIZE - 1] = m_dptable[(i - 1) * m_m1 + (j + V_SIZE - 1)] + costs(m_seq1[i - 1], m_gapChar);
         _b = SSECMP16(t.s, cd.s);
         cd.s = _mm_or_si128(_mm_and_si128(_b, t.s), _mm_andnot_si128(_b, cd.s));
         cb.s = _mm_or_si128(_mm_and_si128(_b, cbd.s), _mm_andnot_si128(_b, cb.s));

         /* REP */
         for (uint32_t k = 0; k < (V_SIZE - 1); k++) {
            t.c[k] = pp.c[k + 1] + costs(m_seq1[i + (V_SIZE - 1 - k) - 1], m_seq2[j + k - 1]);
         }
         t.c[V_SIZE - 1] = m_dptable[(i - 1) * m_m1 + (j + V_SIZE - 1 - 1)] + costs(m_seq1[i - 1], m_seq2[j + V_SIZE - 1 - 1]);
         _b = SSECMP16(t.s, cd.s);
         cd.s = _mm_or_si128(_mm_and_si128(_b, t.s), _mm_andnot_si128(_b, cd.s));
         cb.s = _mm_or_si128(_mm_and_si128(_b, cbr.s), _mm_andnot_si128(_b, cb.s));

         /* write back */
         for (uint32_t k = 0; k < V_SIZE; k++) {
            m_dptable[(i + (V_SIZE - 1 - k)) * m_m1 + j + k] = cd.c[k];
            m_betable[(i + (V_SIZE - 1 - k)) * m_m1 + j + k] = (EdgeType)cb.e[k];
            if (m_dptable[(i + (V_SIZE - 1 - k)) * m_m1 + j + k] > m_smax) {
               m_smax = m_dptable[(i + (V_SIZE - 1 - k)) * m_m1 + j + k];
               m_imax = i + k;
               m_jmax = j;
            }
         }
      }
   }
#endif

#if 1
   /* init last V_SIZE row -------------------------- */
   for (uint32_t i = endi; i < m_n1; i++) {
      for (uint32_t j = 2; j < endj; j++) {
         /* short cut to position of matrix we want to fill */
         dpcurr = &m_dptable[i * m_m1 + j];
         becurr = &m_betable[i * m_m1 + j];

         /* gap in seq1 */
         tmp_costs = m_dptable[i * m_m1 + (j - 1)]
               + costs(m_gapChar, m_seq2[j - 1]);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_INS;
         }
         /* gap in seq2 */
         tmp_costs = m_dptable[(i - 1) * m_m1 + j]
               + costs(m_seq1[i - 1], m_gapChar);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_DEL;
         }
         /* replacement */
         tmp_costs = m_dptable[(i - 1) * m_m1 + (j - 1)]
               + costs(m_seq1[i - 1], m_seq2[j - 1]);
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

   /* init last V_SIZE column-------------------------- */
   for (uint32_t i = 1; i < m_n1; i++) {
      for (uint32_t j = endj; j < m_m1; j++) {
         /* short cut to position of matrix we want to fill */
         dpcurr = &m_dptable[i * m_m1 + j];
         becurr = &m_betable[i * m_m1 + j];

         /* gap in seq1 */
         tmp_costs = m_dptable[i * m_m1 + (j - 1)]
               + costs(m_gapChar, m_seq2[j - 1]);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_INS;
         }
         /* gap in seq2 */
         tmp_costs = m_dptable[(i - 1) * m_m1 + j]
               + costs(m_seq1[i - 1], m_gapChar);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_DEL;
         }
         /* replacement */
         tmp_costs = m_dptable[(i - 1) * m_m1 + (j - 1)]
               + costs(m_seq1[i - 1], m_seq2[j - 1]);
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
#endif

   //printTables();

   return m_dptable[m_imax * m_m1 + m_jmax];
}

/*----------------------------------------------------------------------------*/

template<typename _Costs>
__attribute__((optimize("unroll-loops")))
int32_t SeqAligner::fillTable_sse_32(_Costs costs) {
   int32_t tmp_costs;
   int32_t* __restrict__ dpcurr;
   EdgeType* __restrict__ becurr;

   assert(m_gapChar != m_seq1[0]);
   assert(m_gapChar != m_seq2[0]);
   m_dptable[0] = 0;
   m_betable[0] = EDGETYPE_REP;

   const uint32_t V_SIZE = 4u;


   /* init first row ------------------------- */
   for (uint32_t j = 1; j < m_m1; j++) {
      /* gap in seq2 */
      m_dptable[0 + j] = m_dptable[0 + (j - 1)]
            + costs(m_seq1[j - 1], m_gapChar);
      m_betable[0 + j] = EDGETYPE_DEL;
      /* start of new alignment? */
      m_dptable[0 + j] = std::max(m_dptable[0 + j], 0);
      /* end of local alignment */
      if (m_dptable[0 + j] > m_smax) {
         m_smax = m_dptable[0 + j];
         m_imax = 0;
         m_jmax = j;
      }
   }

   /* init first V_SIZE columns-------------------------- */
   for (uint32_t i = 1; i < m_n1; i++) {
      /* gap in seq1 */
      m_dptable[i * m_m1] = m_dptable[(i - 1) * m_m1]
            + costs(m_gapChar, m_seq2[i - 1]);
      m_betable[i * m_m1] = EDGETYPE_INS;
      /* start of new alignment? */
      m_dptable[i * m_m1] = std::max(m_dptable[i * m_m1], 0);
      if (m_dptable[i * m_m1] > m_smax) {
         m_smax = m_dptable[i * m_m1];
         m_imax = 0;
         m_jmax = i;
      }
      for (uint32_t j = 1; j < std::min(m_m1, V_SIZE) + 2; j++) {
         /* short cut to position of matrix we want to fill */
         dpcurr = &m_dptable[i * m_m1 + j];
         becurr = &m_betable[i * m_m1 + j];

         /* gap in seq1 */
         tmp_costs = m_dptable[i * m_m1 + (j - 1)]
               + costs(m_gapChar, m_seq2[j - 1]);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_INS;
         }
         /* gap in seq2 */
         tmp_costs = m_dptable[(i - 1) * m_m1 + j]
               + costs(m_seq1[i - 1], m_gapChar);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_DEL;
         }
         /* replacement */
         tmp_costs = m_dptable[(i - 1) * m_m1 + (j - 1)]
               + costs(m_seq1[i - 1], m_seq2[j - 1]);
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

   /* use diagonal */
#if 1

   union _TD {
      __m128i s;
      int32_t c[V_SIZE];

   } pp, p, t, cd;

   union _TE {
      __m128i s;
      int32_t e[V_SIZE];

   };
   const _TE cbi = {_mm_set1_epi32(EDGETYPE_INS)};
   const _TE cbd = {_mm_set1_epi32(EDGETYPE_DEL)};
   const _TE cbr = {_mm_set1_epi32(EDGETYPE_REP)};

   _TE cb;

   __m128i _b;
   const __m128i  _s0 = _mm_set1_epi32(0);
#define SSECMP32 _mm_cmpgt_epi32


   //char t[V_SIZE];
   uint32_t endi = m_n1 - ((m_n1 - 1) % V_SIZE);
   //endi = std::max(endi, V_SIZE) - V_SIZE;
   uint32_t endj = m_m1 - (V_SIZE);
   for (uint32_t i = 1; i < endi; i += V_SIZE) {
      for (uint32_t j = 2; j < endj; j += 1) {
         if (j == 2) {
            for (uint32_t k = 0; k < V_SIZE; k++) {
               pp.c[k] = m_dptable[(i + (V_SIZE - 1 - k)) * m_m1 + (j + k - 2)];
               p.c[k] = m_dptable[(i + (V_SIZE - 1 - k)) * m_m1 + (j + k - 1)];
            }
         } else {
            {
               pp = p;
            }
            for (uint32_t k = 0; k < V_SIZE; k++) {
               p.c[k] = m_dptable[(i + (V_SIZE - 1 - k)) * m_m1 + (j + k - 1)];
            }
         }
         cb.s = cbi.s;
         cd.s = _s0;

         /* INS */
         for (uint32_t k = 0; k < (V_SIZE - 0); k++) {
            t.c[k] = p.c[k] + costs(m_gapChar, m_seq2[j + k - 1]);
         }
         _b = SSECMP32(t.s, cd.s);
         cd.s = _mm_or_si128(_mm_and_si128(_b, t.s), _mm_andnot_si128(_b, cd.s));
         cb.s = _mm_or_si128(_mm_and_si128(_b, cbi.s), _mm_andnot_si128(_b, cb.s));


         /* DEL */
         for (uint32_t k = 0; k < (V_SIZE - 1); k++) {
            t.c[k] = p.c[k + 1] + costs(m_seq1[i + (V_SIZE - 1 - k) - 1], m_gapChar);
         }
         t.c[V_SIZE - 1] = m_dptable[(i - 1) * m_m1 + (j + V_SIZE - 1)] + costs(m_seq1[i - 1], m_gapChar);
         _b = SSECMP32(t.s, cd.s);
         cd.s = _mm_or_si128(_mm_and_si128(_b, t.s), _mm_andnot_si128(_b, cd.s));
         cb.s = _mm_or_si128(_mm_and_si128(_b, cbd.s), _mm_andnot_si128(_b, cb.s));

         /* REP */
         for (uint32_t k = 0; k < (V_SIZE - 1); k++) {
            t.c[k] = pp.c[k + 1] + costs(m_seq1[i + (V_SIZE - 1 - k) - 1], m_seq2[j + k - 1]);
         }
         t.c[V_SIZE - 1] = m_dptable[(i - 1) * m_m1 + (j + V_SIZE - 1 - 1)] + costs(m_seq1[i - 1], m_seq2[j + V_SIZE - 1 - 1]);
         _b = SSECMP32(t.s, cd.s);
         cd.s = _mm_or_si128(_mm_and_si128(_b, t.s), _mm_andnot_si128(_b, cd.s));
         cb.s = _mm_or_si128(_mm_and_si128(_b, cbr.s), _mm_andnot_si128(_b, cb.s));

         /* write back */
         for (uint32_t k = 0; k < V_SIZE; k++) {
            m_dptable[(i + (V_SIZE - 1 - k)) * m_m1 + j + k] = cd.c[k];
            m_betable[(i + (V_SIZE - 1 - k)) * m_m1 + j + k] = (EdgeType)cb.e[k];
            if (m_dptable[(i + (V_SIZE - 1 - k)) * m_m1 + j + k] > m_smax) {
               m_smax = m_dptable[(i + (V_SIZE - 1 - k)) * m_m1 + j + k];
               m_imax = i + k;
               m_jmax = j;
            }
         }
      }
   }
#endif

#if 1
   /* init last V_SIZE row -------------------------- */
   for (uint32_t i = endi; i < m_n1; i++) {
      for (uint32_t j = 2; j < endj; j++) {
         /* short cut to position of matrix we want to fill */
         dpcurr = &m_dptable[i * m_m1 + j];
         becurr = &m_betable[i * m_m1 + j];

         /* gap in seq1 */
         tmp_costs = m_dptable[i * m_m1 + (j - 1)]
               + costs(m_gapChar, m_seq2[j - 1]);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_INS;
         }
         /* gap in seq2 */
         tmp_costs = m_dptable[(i - 1) * m_m1 + j]
               + costs(m_seq1[i - 1], m_gapChar);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_DEL;
         }
         /* replacement */
         tmp_costs = m_dptable[(i - 1) * m_m1 + (j - 1)]
               + costs(m_seq1[i - 1], m_seq2[j - 1]);
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

   /* init last V_SIZE column-------------------------- */
   for (uint32_t i = 1; i < m_n1; i++) {
      for (uint32_t j = endj; j < m_m1; j++) {
         /* short cut to position of matrix we want to fill */
         dpcurr = &m_dptable[i * m_m1 + j];
         becurr = &m_betable[i * m_m1 + j];

         /* gap in seq1 */
         tmp_costs = m_dptable[i * m_m1 + (j - 1)]
               + costs(m_gapChar, m_seq2[j - 1]);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_INS;
         }
         /* gap in seq2 */
         tmp_costs = m_dptable[(i - 1) * m_m1 + j]
               + costs(m_seq1[i - 1], m_gapChar);
         if (tmp_costs >= dpcurr[0]) {
            dpcurr[0] = tmp_costs;
            becurr[0] = EDGETYPE_DEL;
         }
         /* replacement */
         tmp_costs = m_dptable[(i - 1) * m_m1 + (j - 1)]
               + costs(m_seq1[i - 1], m_seq2[j - 1]);
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
#endif

   //printTables();

   return m_dptable[m_imax * m_m1 + m_jmax];
}
#endif

#endif /*SEQALIGNER_H_*/
