#include <stdint.h>
#include <vector>
#include <algorithm>
#include <assert.h>


namespace mass {

/*----------------------------------------------------------------------------*/

/* redefine to use algorithm */
class CompU {
public:
   virtual ~CompU() {};
   float operator()(const uint32_t a, const uint32_t b) {return a < b;}
};

/*----------------------------------------------------------------------------*/

template <typename T>
struct MatchSetItem {
   const T* s1;
   const T* s2;
};

template <typename T>
struct Match {
   float score;
   uint32_t numTries;
   std::vector<MatchSetItem<T> > matches;
   Match() {
      score = 0;
      numTries = 0;
   }
};

/*----------------------------------------------------------------------------*/

template <typename T> class Person;

template <typename T>
class Suitor {
public:
   Person<T>* p;
   float score;
   bool proposed;
};

/*----------------------------------------------------------------------------*/

template <typename T>
class Person {
public:

   Person(const T &u);
   template <typename _Cmp>
   void eval(std::vector<Person> &s1, _Cmp &cmp);
   std::vector<Suitor<T> > suitors;
   Person* fiance;
   const T* m_u;
private:
};

template <typename T>
Person<T>::Person(const T &u) {
   m_u = &u;
   fiance = NULL;
}

template <typename T>
static bool cmpBigger(const Suitor<T>& a, const Suitor<T>& b) {
   return a.score > b.score;
}

template <typename T>
template <typename _Cmp>
void Person<T>::eval(std::vector<Person<T> > &s1, _Cmp &cmp) {
   for (uint32_t i = 0; i < s1.size(); ++i) {
      suitors.push_back(Suitor<T>());
      suitors.back().p = &s1[i];
      suitors.back().score = cmp(m_u, s1[i].m_u);
      suitors.back().proposed = false;
   }
   std::sort(suitors.begin(), suitors.end(), cmpBigger<T>);
}

/*----------------------------------------------------------------------------*/

template <typename T>
class Engagemens {
public:
   void free(Person<T> *m, Person<T> *w) {
      typename std::map<std::pair<Person<T> *, Person<T> *>, float >::iterator it;
      it = m_matches.find(std::make_pair(m, w));
      assert(it != m_matches.end());
      m_matches.erase(it);
   }

   float score(Person<T> * m, Person<T> * w) const {
      typename std::map<std::pair<Person<T> *, Person<T> *>, float >::const_iterator it;
      it = m_matches.find(std::make_pair(m, w));
      return it != m_matches.end() ? it->second : 0;
   }

   void create(Person<T> * m, Person<T> * w, float score) {
      m_matches[std::make_pair(m, w)] = score;
   }

   bool exists(Person<T> * m, Person<T> * w) const {
      typename std::map<std::pair<Person<T> *, Person<T> *>, float >::const_iterator it;
      it = m_matches.find(std::make_pair(m, w));
      return it != m_matches.end();
   }

   void getRes(Match<T> &res) {
      res.matches.clear();
      res.score = 0;
      typename std::map<std::pair<Person<T> *, Person<T> *>, float >::const_iterator it;
      for(it = m_matches.begin(); it != m_matches.end(); ++it) {
         res.matches.push_back(MatchSetItem<T>());
         res.matches.back().s1 = it->first.first->m_u;
         res.matches.back().s2 = it->first.second->m_u;
         res.score += it->second;
      }
   }
private:
   typename std::map<std::pair<Person<T> *, Person<T> *>, float > m_matches;
};

/*----------------------------------------------------------------------------*/

template <typename T, typename _Cmp >
void stableMarriageMatching(const std::vector<T> &s1, const std::vector<T> &s2, _Cmp &cmp, Match<T> &res) {
   std::vector<Person<T> > sM;
   std::vector<Person<T> > sW;
   /* Initialize all m ∈ M and w ∈ W to free */
   for (uint32_t i = 0; i < s1.size(); i++) {
      sM.push_back(Person<T>(s1[i]));
   }
   for (uint32_t i = 0; i < s2.size(); i++) {
      sW.push_back(Person<T>(s2[i]));
   }

   /* set preferences of m to w and w to m */
   for (uint32_t i = 0; i < sM.size(); i++) {
      sM[i].eval(sW, cmp);
   }
   for (uint32_t i = 0; i < sW.size(); i++) {
      sW[i].eval(sM, cmp);
   }

   Engagemens<T> engagements;

   uint32_t numMatched = 0;
   /* algorithm */
   while (numMatched < s1.size() && numMatched < s2.size()) {
      for (uint32_t i = 0; i < sM.size(); i++) {
         if (sM[i].fiance != NULL) {
            continue;
         }
         Person<T> *m = &sM[i];
         /* m's highest ranked woman */
         for (uint32_t j = 0; j < m->suitors.size(); ++j) {
            /* to whom he has not yet proposed */
            if (m->suitors[j].proposed == true) {
               continue;
            }
            /* don't try again */
            m->suitors[j].proposed = true;

            Suitor<T> &msw = sM[i].suitors[j];
            Person<T> *w = sM[i].suitors[j].p;
            /* if w is free */
            if (w->fiance == NULL) {
               /* (m, w) become engaged */
               w->fiance = m;
               m->fiance = w;
               engagements.create(m, w, msw.score);
               numMatched++;
            } else {
               assert(engagements.exists(w->fiance, w) && "w has to be engaged");
               /* if w prefers m to m' */
               if (msw.score > engagements.score(w->fiance, w)) {
                  engagements.free(w->fiance, w);
                  w->fiance->fiance = NULL;
                  w->fiance = m;
                  m->fiance = w;
                  engagements.create(m, w, msw.score);
               } else {
                  /* (m', w) remain engaged */
               }
            }
            break;
         }
      }
   }
   /* res is a score and a list of pairwise matches */
   engagements.getRes(res);
}

/*----------------------------------------------------------------------------*/
#include <limits.h>
template <typename T, typename _Cmp>
void maximumWeightedMatching(const std::vector<T> &s1, const std::vector<T> &s2, _Cmp &cmp, mass::Match<T> &res) {
   const bool switched = s1.size() > s2.size();
   const std::vector<T> &i1 = switched == false ? s1 : s2;
   const std::vector<T> &i2 = switched == false ? s2 : s1;
#define MATCH_UNASSIGNED UINT_MAX
   class Match {
   public:

      /* try all mappings from i1 -> i2 */
      void matchRec(uint32_t depth, std::vector<uint32_t> &m, float score, _Cmp &cmp) {
         /* assigned all i1 */
         if (depth >= i1.size()) {
            if (score > scoreMax) {
               resMax = m;
               scoreMax = score;
            }
            return;
         }
         for (uint32_t i = 0; i < m.size(); i++) {
            if (m[i] == MATCH_UNASSIGNED) {
               m[i] = i;
               float scoreNew = score + cmp(&i1[depth], &i2[m[i]]);
               matchRec(depth + 1, m, scoreNew, cmp);
               m[i] = MATCH_UNASSIGNED;
            }
         }
      }

      Match(const std::vector<T>& s1, const std::vector<T> &s2, _Cmp &cmp) : i1(s1), i2(s2) {
         scoreMax = 0;
         std::vector<uint32_t> m;
         m.resize(s1.size(), MATCH_UNASSIGNED);
         matchRec(0, m, 0.f, cmp);
      }

      std::vector<uint32_t> resMax;
      float scoreMax;
      const std::vector<T> &i1;
      const std::vector<T> &i2;
   };

   Match a(i1, i2, cmp);

   res.score = a.scoreMax;
   for (uint32_t i = 0; i < a.resMax.size(); ++i) {
      if (a.resMax[i] != MATCH_UNASSIGNED) {
         res.matches.push_back(mass::MatchSetItem<T>());
         if (switched) {
            res.matches.back().s2 = &a.i1[a.resMax[i]];
            res.matches.back().s1 = &a.i2[i];
         } else {
            res.matches.back().s1 = &a.i1[a.resMax[i]];
            res.matches.back().s2 = &a.i2[i];
         }
      }
   }
#undef MATCH_UNASSIGNED
}


} /* namespace sma */


