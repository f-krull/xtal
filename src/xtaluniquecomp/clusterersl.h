

#ifndef CLUSTERERSL_H_
#define CLUSTERERSL_H_

#include "clusterer.h"
#include "dendrogram.h"

template<class T>
class ClustererSl: public Clusterer<T> {
public:
   ClustererSl(DistanceMatrix<T> *dm, T maxDist);
   bool linkNext(unsigned int *i, unsigned int *j, T *dist);
   static unsigned int clusterAll(DistanceMatrix<T> *dm,
         std::vector<std::vector<uint32_t> > *clusters, T maxDist,
         uint32_t maxSteps = 0, const std::vector<std::string> *eNames = NULL);
private:
   /* union-find data struct */
   class UfNode {
   public:
      UfNode *head;
      UfNode *next;
      unsigned int clusterId;
   };
   /* edge list entry of MST */
   class MstEdge {
   public:
      T dist;
      unsigned int node1;
      unsigned int node2;
      MstEdge(unsigned int node1, unsigned int node2, T dist);
   };
   /* union-find nodes */
   std::vector<UfNode> ufNodes;
   /* MST edges */
   std::vector<MstEdge> edgeList;
   unsigned int currentEdge;
   unsigned int numClusters;

   /* used for sort(); smaller distances first */
   static bool compareEdges(const MstEdge &a, const MstEdge &b);
   /* union two UF sets */
   void unionTrees(UfNode *node1, UfNode *node2, unsigned int newId);

};

#include <iostream>
#include <iomanip>
#include <cmath>
#include <climits>

/*----------------------------------------------------------------------------*/

template <class T>
bool ClustererSl<T>::compareEdges(const MstEdge &a, const MstEdge &b) {
   return a.dist < b.dist;
}

/*----------------------------------------------------------------------------*/

template <class T>
ClustererSl<T>::MstEdge::MstEdge(unsigned int n1, unsigned int n2, T dist) {
   this->node1 = n1;
   this->node2 = n2;
   this->dist = dist;
}

/*----------------------------------------------------------------------------*/

template <class T>
ClustererSl<T>::ClustererSl(DistanceMatrix<T>* dm, T maxDist) {
   currentEdge = 0;
   /* create edge list */
   for (unsigned int i = 0; i < dm->getNumElements(); i++) {
      for (unsigned int j = i + 1; j < dm->getNumElements(); j++) {
         if ((maxDist != 0) && (dm->get(i, j) > maxDist)) {
            continue;
         }
         edgeList.push_back(MstEdge(i, j, dm->get(i, j)));
      }
   }
   /* sort edges by distance */
   sort(edgeList.begin(), edgeList.end(), compareEdges);
   /* create union-find sets */
   ufNodes.resize(dm->getNumElements());
   for (unsigned int i = 0; i < ufNodes.size(); i++) {
      ufNodes[i].head = &ufNodes[i];
      ufNodes[i].next = NULL;
      ufNodes[i].clusterId = i;
   }
   numClusters = dm->getNumElements();
}

/*----------------------------------------------------------------------------*/

template <class T>
void ClustererSl<T>::unionTrees(UfNode *node1, UfNode *node2, unsigned int id) {
   /* reach end of first set */
   while (node1->next != NULL) {
      node1 = node1->next;
   }
   /* connect sets */
   node2 = node2->head;
   node1->next = node2;
   node2->head = node1->head;
   /* heads of second set are head of first set */
   while (node2->next != NULL) {
      node2 = node2->next;
      node2->head = node1->head;
   }
   /* set id of head */
   node1->head->clusterId = id;
}

/*----------------------------------------------------------------------------*/

template <class T>
bool ClustererSl<T>::linkNext(unsigned int *i, unsigned int *j, T *dist) {
   unsigned int node1, node2;

   /* find smallest edge */
   while (currentEdge < edgeList.size()) {
      node1 = edgeList[currentEdge].node1;
      node2 = edgeList[currentEdge].node2;
      /* node has to connect two unconnected trees */
      if (ufNodes[node1].head != ufNodes[node2].head) {
         /*edge found; store result */
         (*dist) = edgeList[currentEdge].dist;
         /* get cluster Ids */
         (*i) = ufNodes[node1].head->clusterId;
         (*j) = ufNodes[node2].head->clusterId;
         /* connect trees and set cluster IDs  */
         unionTrees(&ufNodes[node1], &ufNodes[node2], numClusters);
         /* edge was processed */
         currentEdge++;
         numClusters++;
         return true;
      }
      currentEdge++;
   }
   /* no edge found */
   (*i) = 0;
   (*j) = 0;
   (*dist) = (T)0;
   return false;
}

/*----------------------------------------------------------------------------*/

template <class T>
unsigned int ClustererSl<T>::clusterAll(DistanceMatrix<T> *dm,
      std::vector<std::vector<uint32_t> > *clusters, T maxDist,
      uint32_t maxSteps, const std::vector<std::string> *eNames) {
   Clusterer<T>* cl;
      uint32_t numPairs, i, j, n;
      T dist;

      n = dm->getNumElements();
      Log::dbg("cl.clustering of %u elements", n);
   #if 0
      UnionFindSet ufs(n);
   #endif
      /* map tree nodes to ids (tracks ids of clustered elements) */
      std::vector<uint32_t> dict(n * 2);
      for (uint32_t i = 0; i < dict.size(); i++) {
         dict[i] = i;
      }

      if (maxSteps == 0) {
         maxSteps = UINT_MAX;
      }
      /* do clustering */
      numPairs = 0;
      cl = new ClustererSl(dm, maxDist);

      Dendrogram dgr;
      dgr.init(*eNames);

      while ((cl->linkNext(&i, &j, &dist)) && (numPairs < maxSteps)) {
   #if 0
            /* still clustering or just generating output? */
            if (ufs.join(dict[i], dict[j]) == false) {
               Log::err("cl.clustering error");
               exit(1);
            }
   #endif
         dgr.pair(i, j, dist);
         dict[numPairs + n] = dict[i];
         numPairs++;

      }
      delete cl;
      (*clusters) = dgr.getClusters();
   #if 0
      //ufs.getSets(clusters);
   #endif
      Log::dbg("cl.single-linkageclustering done: %u links, %u sets",
            numPairs, clusters->size());
      return clusters->size();
}


#endif /* CLUSTERERSL_H_ */
