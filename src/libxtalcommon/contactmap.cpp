#include "contactmap.h"
#include <string.h>
#include "../libxtaldata/aamap.h"
#include "../libxtalutil/log.h"
#include <stdlib.h>

ContactMap::ContactMap() {
   m_map.resize(AaMap::size * AaMap::size);
   m_mapt.resize(AaMap::size * AaMap::size);
   m_numContacts = 0;
}

ContactMap::~ContactMap() {
}

void ContactMap::build(const Chains &c1, const Chains &c2, float cutoff) {
   uint32_t numContacts = 0;
   for (uint32_t i = 0; i < c1.size(); i++) {
      for (uint32_t j = 0; j < c2.size(); j++) {
         numContacts += add(c1[i]->resis(), c2[j]->resis(), cutoff);
      }
   }
   m_numContacts = numContacts;
   Log::dbg("numContacts: %u", numContacts);
}

uint32_t ContactMap::add(const Residues &r1, const Residues &r2, float cutoff) {
   uint32_t numContacts = 0;
   for (uint32_t i = 0; i < r1.size(); i++) {
      uint16_t id1 = AaMap::getId(r1[i]->getName());
      if (id1 == 0) {
         continue;
      }
      for (uint32_t j = 0; j < r2.size(); j++) {
         uint16_t id2 = AaMap::getId(r2[j]->getName());
         if (id2 == 0) {
            continue;
         }
         if (Residue::hasContact_fast(r1[i], r2[j], cutoff)) {
            m_map[id1 * AaMap::size + id2]++;
            m_mapt[id2 * AaMap::size + id1]++;
            numContacts++;
         }
      }
   }
   return numContacts;
}

float ContactMap::dist(const ContactMap &c2) const {
   const ContactMap &c1 = (*this);

   const ContactMap &b = c1.m_numContacts > c2.m_numContacts ? c1 : c2;
   const ContactMap &s = c1.m_numContacts > c2.m_numContacts ? c2 : c1;

   float dist1 = 0;
   float dist2 = 0;


   {
      int32_t count = 0;
      int32_t diff = 0;
      for (uint32_t i = 1; i < (uint32_t)AaMap::size; i++) {
         for (uint32_t j = 1; j < (uint32_t)AaMap::size; j++) {
            uint32_t mb = b.m_map[i * AaMap::size + j];
            uint32_t ms = s.m_map[i * AaMap::size + j];
            if (ms < mb) {
               diff += labs(ms - mb);
            }
            count += mb;
            count += ms;
         }
      }
      dist1 = (float(diff) / count);
   }
   {
     int32_t count = 0;
     int32_t diff = 0;
     for (uint32_t i = 1; i < (uint32_t)AaMap::size; i++) {
        for (uint32_t j = 1; j < (uint32_t)AaMap::size; j++) {
           uint32_t mb = b.m_mapt[i * AaMap::size + j];
           uint32_t ms = s.m_map[i * AaMap::size + j];
           if (ms < mb) {
              diff += labs(ms - mb);
           }
           count += mb;
           count += ms;
        }
     }
     dist2 = (float(diff) / count);
  }

   return std::min(dist1, dist2);
}

#define POW_2(x) ((x) * (x))
#define POW_3(x) ((x) * (x) * (x))

const uint32_t ContactMap3::NUM_MAP_ROWS = POW_2(AaMap::size);

ContactMap3::ContactMap3() {
   uint32_t mapSize = POW_2(NUM_MAP_ROWS);
   m_map.resize(mapSize);
   m_mapt.resize(mapSize);
}

uint32_t ContactMap3::add(const Residues &r1, const Residues &r2, float cutoff) {
   uint32_t numContacts = 0;
   for (uint32_t i = 1; i < r1.size(); i++) {
      uint32_t id1_0 = AaMap::getId(r1[i-1]->getName());
      uint32_t id1_1 = AaMap::getId(r1[i+0]->getName());
      if (id1_1 == 0) {
         continue;
      }
      uint32_t id1 = 0;
      id1 += id1_1 * (AaMap::size);
      id1 += id1_0 * (1);
      for (uint32_t j = 1; j < r2.size(); j++) {
         uint32_t id2_0 = AaMap::getId(r2[j-1]->getName());
         uint32_t id2_1 = AaMap::getId(r2[j+0]->getName());
         if (id2_1 == 0) {
            continue;
         }
         uint32_t id2 = 0;
         id2 += id2_1 * (AaMap::size);
         id2 += id2_0 * (1);
         if (Residue::hasContact_fast(r1[i], r2[j], cutoff)) {
            m_map[id1 * NUM_MAP_ROWS + id2]++;
            m_mapt[id2 * NUM_MAP_ROWS + id1]++;
            numContacts++;
         }
      }
   }
   return numContacts;
}

float ContactMap3::dist(const ContactMap3 &c2) const {
   const ContactMap3 &c1 = (*this);
   const ContactMap3 &b = c1.m_numContacts > c2.m_numContacts ? c1 : c2;
   const ContactMap3 &s = c1.m_numContacts > c2.m_numContacts ? c2 : c1;

   float dist1 = 0;
   float dist2 = 0;

   const uint32_t mapSize = POW_2(NUM_MAP_ROWS);

   {
      int32_t count = 0;
      int32_t diff = 0;
      for (uint32_t i = 0; i < mapSize; i++) {
         uint32_t mb = b.m_map[i];
         uint32_t ms = s.m_map[i];
         if (ms < mb) {
            diff += labs(ms - mb);
         }
         count += mb;
         count += ms;
      }
      dist1 = (float(diff) / count);
   }
   {
      int32_t count = 0;
      int32_t diff = 0;
      for (uint32_t i = 0; i < b.m_map.size(); i++) {
         uint32_t mb = b.m_mapt[i];
         uint32_t ms = s.m_map[i];
         if (ms < mb) {
            diff += labs(ms - mb);
         }
         count += mb;
         count += ms;
      }
      dist2 = (float(diff) / count);
   }

   return std::min(dist1, dist2);
}

