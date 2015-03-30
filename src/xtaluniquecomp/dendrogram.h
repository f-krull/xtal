#ifndef DENDROGRAM_H_
#define DENDROGRAM_H_

#include <string>
#include <vector>
#include <stdint.h>

/*----------------------------------------------------------------------------*/

class DeNode;

/*----------------------------------------------------------------------------*/

class Dendrogram {
public:

   Dendrogram();

   ~Dendrogram();

   void init(const std::vector<std::string> &nodeLables);

   bool pair(unsigned int i, unsigned int j, float dist);


   std::vector<std::vector<uint32_t> > getClusters() const;

private:
   std::vector<std::string> m_nodeLables;
   std::vector<uint32_t> m_dict;
   std::vector<DeNode*> m_nodes;
   uint32_t m_n;
   uint32_t m_numPairs;
};

#endif /* DENDROGRAM_H_ */
