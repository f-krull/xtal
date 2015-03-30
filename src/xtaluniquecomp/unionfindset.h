#ifndef UNIONFINDSET_H_
#define UNIONFINDSET_H_

#include <vector>
#include <stdint.h>

/*----------------------------------------------------------------------------*/

class UnionFindSet {
private:

   class UfNode {
   public:
      UfNode *head;
      UfNode *next;
      uint32_t id;

      UfNode();
   };

   std::vector<UfNode> elements;
   uint32_t numElements;

protected:
public:
   UnionFindSet(uint32_t setSize); /* indices are form 0 to setSize - 1 */
   bool join(uint32_t i, uint32_t j);
   bool sameSet(uint32_t i, uint32_t j);
   void getSets(std::vector<std::vector<uint32_t> >* sets);
};

#endif /*UNIONFINDSET_H_*/
