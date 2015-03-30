
#include <algorithm>

/*----------------------------------------------------------------------------*/

static std::string chainsToStr(const Chains chains, bool sort = true) {
   std::string res;
   for (uint32_t i = 0; i < chains.size(); i++) {
      res.push_back(chains[i]->name());
   }
   if (sort == true) {
      std::sort(res.begin(), res.end());
   }
   return res;
}


