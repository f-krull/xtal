#include "dendrogram.h"
#include "../libxtalutil/log.h"
#include <queue>
#include <map>
#include <assert.h>
#include <algorithm>
#include <stdlib.h>

/*----------------------------------------------------------------------------*/

class DeNode {
public:

   DeNode(int32_t id, const std::string &label);

   bool isRoot() const;
   bool isLeaf() const;
   int32_t id() const;

   static DeNode * merge(int32_t id, DeNode *l, DeNode *r, float dist);

   static void getSteps(uint32_t clusterId, DeNode *root,
         const std::map<int32_t, int32_t> &nodeMap);

   static void getNodes(uint32_t clusterId, DeNode *root,
         const std::map<int32_t, int32_t> &nodeMap);

   static uint32_t getNumLeafs(const DeNode *root);

   static void printOrder(uint32_t clusterId, const DeNode *root,
         const std::map<int32_t, int32_t> &nodeMap);

   static std::vector<int32_t> getNodeIds(const DeNode *root);

private:
   DeNode *m_p;
   DeNode *m_l;
   DeNode *m_r;
   float m_dist;
   bool m_isLeaf;
   std::string m_label;
   int32_t m_id;
   uint32_t m_depth;

   class PrioIdRevComp {
   public:
      bool operator()(const DeNode* lhs, const DeNode* rhs) const {
         return lhs->id() > rhs->id();
      }
   };

   DeNode(int32_t id, DeNode *l, DeNode *r, float dist);

   static int32_t mapId(const std::map<int32_t, int32_t> &nodeMap
         , int32_t key);

   void printLeaf(uint32_t clusterId,
         const std::map<int32_t, int32_t> &nodeMap) const;

   void printOrder(uint32_t clusterId,
         const std::map<int32_t, int32_t> &nodeMap) const;

   void printStep(uint32_t clusterId,
         const std::map<int32_t, int32_t> &nodeMap) const;

   void traverseSteps(
         std::priority_queue<DeNode*, std::vector<DeNode*>, PrioIdRevComp> &list);
   void traverseNodes(
         std::priority_queue<DeNode*, std::vector<DeNode*>, PrioIdRevComp> &list);

   uint32_t getNumLeafs() const;

   void getNodeIds(std::vector<int32_t> &nodeIds) const;

};

/*----------------------------------------------------------------------------*/

DeNode::DeNode(int32_t id, const std::string &label) {
   m_label = label;
   m_isLeaf = true;
   m_id = id;
   m_p = NULL;
   m_r = NULL;
   m_l = NULL;
   m_dist = 0;
   m_depth = 0;
}

/*----------------------------------------------------------------------------*/

bool DeNode::isRoot() const {
   return m_p == NULL;
}

/*----------------------------------------------------------------------------*/

bool DeNode::isLeaf() const {
   return m_isLeaf;
}

/*----------------------------------------------------------------------------*/

int32_t DeNode::id() const {
   return m_id;
}

/*----------------------------------------------------------------------------*/

DeNode * DeNode::merge(int32_t id, DeNode *l, DeNode *r, float dist) {
   DeNode *p = NULL;
   {
      /* which one is bigger? */
      DeNode *huge = NULL;
      DeNode *tiny = NULL;
      if (l->m_depth > r->m_depth) {
         huge = l;
         tiny = r;
      } else {
         huge = r;
         tiny = l;
      }
      /* does he want left or right side? */
      if (huge->m_dist > huge->m_dist) {
         p = new DeNode(id, huge, tiny, dist);
      } else {
         p = new DeNode(id, tiny, huge, dist);
      }
   }
   p->m_l->m_p = p;
   p->m_r->m_p = p;
   return p;
}

/*----------------------------------------------------------------------------*/

void DeNode::getSteps(uint32_t clusterId, DeNode *root,
      const std::map<int32_t, int32_t> &nodeMap) {
   std::priority_queue<DeNode*, std::vector<DeNode*>, PrioIdRevComp> next;
   root->traverseSteps(next);
   while (next.empty() == false) {
      next.top()->printStep(clusterId, nodeMap);
      next.pop();
   }
}

/*----------------------------------------------------------------------------*/

void DeNode::getNodes(uint32_t clusterId, DeNode *root,
      const std::map<int32_t, int32_t> &nodeMap) {
   std::priority_queue<DeNode*, std::vector<DeNode*>, PrioIdRevComp> next;
   root->traverseNodes(next);
   while (next.empty() == false) {
      next.top()->printLeaf(clusterId, nodeMap);
      next.pop();
   }
}

/*----------------------------------------------------------------------------*/

uint32_t DeNode::getNumLeafs(const DeNode *root) {
   return root->getNumLeafs();
}

/*----------------------------------------------------------------------------*/

void DeNode::printOrder(uint32_t clusterId, const DeNode *root,
      const std::map<int32_t, int32_t> &nodeMap) {
   root->printOrder(clusterId, nodeMap);
}

/*----------------------------------------------------------------------------*/

std::vector<int32_t> DeNode::getNodeIds(const DeNode *root) {
   std::vector<int32_t> nodeIds;
   root->getNodeIds(nodeIds);
   std::sort(nodeIds.begin(), nodeIds.end());
   return nodeIds;
}

/*----------------------------------------------------------------------------*/

DeNode::DeNode(int32_t id, DeNode *l, DeNode *r, float dist) {
   m_r = r;
   m_l = l;
   m_isLeaf = false;
   m_p = NULL;
   m_dist = dist;
   m_id = id;
   m_label = l->m_label;
   m_depth = std::max(m_r->m_depth, m_l->m_depth) + 1;
}

/*----------------------------------------------------------------------------*/

int32_t DeNode::mapId(const std::map<int32_t, int32_t> &nodeMap , int32_t key) {
   std::map<int32_t, int32_t>::const_iterator it = nodeMap.find(key);
   assert(it != nodeMap.end() && "all node ids should be in the dictionary");
   return it->second;
}

/*----------------------------------------------------------------------------*/

void DeNode::printLeaf(uint32_t clusterId,
      const std::map<int32_t, int32_t> &nodeMap) const {
   const int32_t id = mapId(nodeMap, m_id);
   Log::inf("cl.node%u \"%s\" %d", clusterId, m_label.c_str(), abs(id));
}

/*----------------------------------------------------------------------------*/

void DeNode::printOrder(uint32_t clusterId,
      const std::map<int32_t, int32_t> &nodeMap) const {
   if (isLeaf() == true) {
      const int32_t id = mapId(nodeMap, m_id);
      Log::inf("cl.order%u %d \"%s\"", clusterId, abs(id), m_label.c_str());
      return;
   }
   m_l->printOrder(clusterId, nodeMap);
   m_r->printOrder(clusterId, nodeMap);
}

/*----------------------------------------------------------------------------*/

void DeNode::printStep(uint32_t clusterId,
      const std::map<int32_t, int32_t> &nodeMap) const {
   const int32_t lid = mapId(nodeMap, m_l->m_id);
   const int32_t rid = mapId(nodeMap, m_r->m_id);
   Log::inf("cl.step%u %5d %5d %f  \"%s\" \"%s\"", clusterId, lid, rid, m_dist,
         m_l->m_label.c_str(), m_r->m_label.c_str());
}

/*----------------------------------------------------------------------------*/

void DeNode::traverseSteps(
      std::priority_queue<DeNode*, std::vector<DeNode*>, PrioIdRevComp> &list) {
   if (this->isLeaf() == false) {
      list.push(this);
      m_l->traverseSteps(list);
      m_r->traverseSteps(list);
   }
}

/*----------------------------------------------------------------------------*/

void DeNode::traverseNodes(
      std::priority_queue<DeNode*, std::vector<DeNode*>, PrioIdRevComp> &list) {
   if (this->isLeaf() == true) {
      list.push(this);
   } else {
      m_l->traverseNodes(list);
      m_r->traverseNodes(list);
   }
}

/*----------------------------------------------------------------------------*/

uint32_t DeNode::getNumLeafs() const {
   if (this->isLeaf()) {
      return 1;
   }
   return m_l->getNumLeafs() + m_r->getNumLeafs();
}

/*----------------------------------------------------------------------------*/

void DeNode::getNodeIds(std::vector<int32_t> &nodeIds) const {
   if (this->isLeaf() == false) {
      m_l->getNodeIds(nodeIds);
      m_r->getNodeIds(nodeIds);
   }
   nodeIds.push_back(this->m_id);
}

/*----------------------------------------------------------------------------*/

Dendrogram::Dendrogram() {
   m_n = 0;
}

/*----------------------------------------------------------------------------*/

Dendrogram::~Dendrogram() {
   for (uint32_t i = 0; i < m_nodes.size(); i++) {
      delete m_nodes[i];
   }
}

/*----------------------------------------------------------------------------*/

void Dendrogram::init(const std::vector<std::string> &nodeLables) {
   m_n = nodeLables.size();
   m_dict.resize(m_n * 2);
   for (uint32_t i = 0; i < m_dict.size(); i++) {
      m_dict[i] = i;
   }
   m_nodes.assign(m_n * 2, NULL);
   for (uint32_t i = 0; i < nodeLables.size(); i++) {
      m_nodes[i] = new DeNode(int32_t(i + 1) * -1, nodeLables[i]);
   }
   m_numPairs = 0;

}

/*----------------------------------------------------------------------------*/

bool Dendrogram::pair(unsigned int i, unsigned int j, float dist) {
   m_dict[m_numPairs + m_n] = m_dict[i];
   m_nodes[m_numPairs + m_n] = DeNode::merge(m_numPairs + 1, m_nodes[i],
         m_nodes[j], dist);
   m_numPairs++;

   return true;

}

/*----------------------------------------------------------------------------*/

std::vector<std::vector<uint32_t> > Dendrogram::getClusters() const {
   std::vector<std::vector<uint32_t> > clusters;

   uint32_t clusterId = 0;
   for (uint32_t i = 0; i < m_nodes.size(); i++) {
      if (m_nodes[i] != NULL && m_nodes[i]->isRoot() == true) {
         /* convert from global nodeids to local nodeids */
         std::map<int32_t, int32_t> nodeMap;
         {
            const std::vector<int32_t> nodeIds = DeNode::getNodeIds(m_nodes[i]);
            const uint32_t clSize = DeNode::getNumLeafs(m_nodes[i]);
            Log::inf("cl.cluster %u size %u", clusterId, clSize);
            clusters.push_back(std::vector<uint32_t>());
            for (uint32_t j = 0; j < nodeIds.size(); j++) {

               /* is leaf? */
               if (j < clSize) {
                  nodeMap[nodeIds[j]] = (j + 1) * -1;
                  /* convert from local nodeids to global nodeids */
                  clusters.back().push_back(abs(nodeIds[j]) - 1);
               } else {
                  nodeMap[nodeIds[j]] = j - clSize + 1;
               }
            }
         }
         DeNode::getNodes(clusterId, m_nodes[i], nodeMap);
         DeNode::printOrder(clusterId, m_nodes[i], nodeMap);
         DeNode::getSteps(clusterId, m_nodes[i], nodeMap);
         clusterId++;
      }
   }
   return clusters;
}
