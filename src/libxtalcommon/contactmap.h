#ifndef CONTACTMAP_H_
#define CONTACTMAP_H_

#include "../libxtaldata/protein.h"
#include <stdint.h>

class ContactMap {
public:
   ContactMap();
   virtual ~ContactMap();

   virtual void build(const Chains &c1, const Chains &c2, float cutoff = 6.5f);

   virtual float dist(const ContactMap &c2) const;

protected:
   virtual uint32_t add(const Residues &c1, const Residues &c2, float cutoff);
   std::vector<uint16_t> m_map;
   std::vector<uint16_t> m_mapt; /* transposed */
   uint32_t m_numContacts;
};



class ContactMap3 : public ContactMap {
public:

   ContactMap3();
   virtual ~ContactMap3(){};

   virtual float dist(const ContactMap3 &c2) const;
protected:
   virtual uint32_t add(const Residues &c1, const Residues &c2, float cutoff);

private:
   static const uint32_t NUM_MAP_ROWS;
};


#endif /* CONTACTMAP_H_ */
