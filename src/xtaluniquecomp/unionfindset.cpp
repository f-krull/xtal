#include "./unionfindset.h"
#include <map>
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace std;

/*----------------------------------------------------------------------------*/

UnionFindSet::UfNode::UfNode() {
   head = NULL;
   next = NULL;
   id = 0;
}

/*----------------------------------------------------------------------------*/

UnionFindSet::UnionFindSet(uint32_t setSize) {
   numElements = setSize;

   elements.resize(numElements);
   for (uint32_t i = 0; i < numElements; i++) {
      elements[i].id = i;
      elements[i].head = &elements[i];
   }
}

/*----------------------------------------------------------------------------*/

bool UnionFindSet::sameSet(uint32_t i, uint32_t j) {
   return (elements[i].head == elements[j].head);
}

/*----------------------------------------------------------------------------*/

bool UnionFindSet::join(uint32_t i, uint32_t j) {
   UfNode *tmp;

   /* already same set? */
   if (sameSet(i , j) == true) {
      return false;
   }
   /* go to end of "i"s set */
   tmp = &elements[i];
   while (tmp->next != NULL) {
      tmp = tmp->next;
   }
   /* connect "i"s end with "j"s head */
   tmp->next = elements[j].head;
   /* change head of "j"s set to "i"s head */
   tmp = elements[j].head;
   while (tmp != NULL) {
      tmp->head = elements[i].head;
      tmp = tmp->next;
   }
   return true;
}

/*----------------------------------------------------------------------------*/

static bool compvectors(vector<uint32_t> v1, vector<uint32_t> v2) {
   return v1.size() > v2.size();
}

/*----------------------------------------------------------------------------*/

void UnionFindSet::getSets(vector<vector<uint32_t> >* sets) {
   map<UfNode*, bool> heads;
   map<UfNode*, bool>::iterator it;
   UfNode* tmp;

   /* set empty? */
   assert(sets->size() == 0);
   for (uint32_t i = 0; i < elements.size(); i++) {
      it = heads.find(elements[i].head);
      /* new one? */
      if (it == heads.end()) {
         heads[elements[i].head] = true;
      }
   }
   uint32_t numheads = 0;
   for (it = heads.begin(); it != heads.end(); it++) {
      //cout << (*it).first << endl;
      numheads++;
      /* fill list to vector */
      sets->push_back(vector<uint32_t>());
      tmp = (*it).first;
      while (tmp != NULL) {
         sets->back().push_back(tmp->id);
         tmp = tmp->next;
      }
   }
   sort(sets->begin(), sets->end(), compvectors);
}
