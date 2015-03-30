#include "xtaluniquecomp.h"
#include "../libxtalcommon/intfdescriptor.h"
#include "../libxtalutil/log.h"
#include "../libxtalutil/command.h"
#include "../libxtalutil/common.h"
#include "../libxtaldata/matrix.cpp"
#include "distmatrix.cpp"
#include "clusterersl.h"
#include "../xtalcompunbound/configcompunbound.h"
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <set>

/*----------------------------------------------------------------------------*/

class XtalUniqueCompPriv {
public:
private:
};

/*----------------------------------------------------------------------------*/

XtalUniqueComp::XtalUniqueComp() {
   m = new XtalUniqueCompPriv();
}

/*----------------------------------------------------------------------------*/

XtalUniqueComp::~XtalUniqueComp() {
   delete m;
}

/*----------------------------------------------------------------------------*/

const char* XtalUniqueComp::getName() const {
   return "XtalUniqueComp";
}

/*----------------------------------------------------------------------------*/
#if 0
bool XtalUniqueComp::read(const std::string& filename, UqEntry &entry) {
   bool ret = true;
   std::ifstream infile;
   std::string line;

   entry.intname = common::removeExtension(common::removePath(filename));
   infile.open(filename.c_str());
   if (infile == false) {
      Log::err("cannot open file \"%s\"", filename.c_str());
      return false;
   }

   while (getline(infile, line)) {
      /* extract interface info */
      UqEntry::parse(line, entry);
   }
   infile.close();

   return ret;
}
#endif

/*----------------------------------------------------------------------------*/
#include <float.h>
static std::vector<float> getMedoidDist(const std::vector<uint32_t> &nodes,
      const DistanceMatrix<float> *dm) {
   struct {
      uint32_t index;
      float dist;
   } best;

   assert(nodes.size() > 0);

   /* get medoid */
   best.dist = FLT_MAX;
   best.index = 0;
   for (uint32_t i = 0; i < nodes.size(); i++) {
      float avgdist = 0;
      for (uint32_t j = 0; j < nodes.size(); j++) {
         avgdist += dm->get(nodes[i], nodes[j]);
      }
      if (avgdist < best.dist) {
         best.dist = avgdist;
         best.index = i;
      }
   }
   std::vector<float> dists;
   dists.reserve(nodes.size());

   /* for all nodes compute distance to medoid */
   for (uint32_t i = 0; i < nodes.size(); i++) {
      dists.push_back(dm->get(nodes[best.index], nodes[i]));
   }
   return dists;
}

/*----------------------------------------------------------------------------*/

int XtalUniqueComp::start() {

   if (cmd.getNumArgs() < 2) {
      Log::err("usage: %s [PDBPATH] [INTF_LIST]", getName());
      Log::err("   clusters interfaces");
      return 1;
   }

   std::string sPath = cmd.getArgStr(0);

   /* get list of interfaces */
   std::vector<std::vector<std::string> > complist;
   if (cmd.getNumArgs() == 2) {
      std::string fnIn = cmd.getArgStr(1);
      /* read input table */
      common::readTableColumns(fnIn, 3, complist);
   } else if (cmd.getNumArgs() == 7) {
      /* put 2 args into table */
      complist.push_back(std::vector<std::string>());
      complist.back().push_back(cmd.getArgStr(1));
      complist.back().push_back(cmd.getArgStr(2));
      complist.back().push_back(cmd.getArgStr(3));
      complist.push_back(std::vector<std::string>());
      complist.back().push_back(cmd.getArgStr(4));
      complist.back().push_back(cmd.getArgStr(5));
      complist.back().push_back(cmd.getArgStr(6));
   }

   if (complist.empty() == true) {
      Log::err("err %d", __LINE__);
      return 1;
   }

   /* read interfaces */
   std::vector<UqEntry> entList;
   for (uint32_t i = 0; i < complist.size(); i++) {
      UqEntry ent;
      const std::string sPB = complist[i][0];
      const std::string sB1 = complist[i][1];
      const std::string sB2 = complist[i][2];
      const std::string fn = sPath + "/" + sPB + ".pdb";

      if (ent.read(fn, sB1, sB2, *ConfigCompUnbound().cmIntfDefinition) == false) {
         Log::dbg("error while parsing '%s'", fn.c_str());
         break;
      }
      entList.push_back(ent);
      if ((i + 1) % 100 == 0) {
         Log::inf("read %u entries", entList.size());
      }
   }


   Log::inf("read %u entries", entList.size());

   /* get distance matrix */
   DistanceMatrix<float> *dm = NULL;
#if 0
   Log::inf("loading distance matrix", entList.size());
   dm = new DistanceMatrix<float>(entList.size());
   if (dm->load("/tmp/dm_120225_19_fin") == false) {
      Log::inf("loading failed... generating new distance matrix",
            entList.size());
      delete dm;
#else
      const bool showAlignment = entList.size() < 3;
   dm = DistanceMatrixFactory<float>::getFilled(&entList,
         DistUqentryCmp(showAlignment));
#endif

   /* cluster elements */
   std::vector<std::vector<uint32_t> > clusters;
   std::vector<std::string> eNames;
   for (uint32_t i = 0; i < entList.size(); i++) {
      eNames.push_back(entList[i].intname);
   }
   Log::inf("clustering %lu elements", dm->getNumElements());
   ClustererSl<float>::clusterAll(dm, &clusters, ConfigCompUnbound().uqClMaxDist, 0, &eNames);

   /* get representative of each cluster */
   Log::inf("sorting %lu clusters", clusters.size());
   for (uint32_t i = 0; i < clusters.size(); i++) {
      printf("cl cluster %u   size: %lu\n", i, clusters[i].size());
      std::vector<UqEntry> entCluster;
      entCluster.reserve(clusters[i].size());
      for (uint32_t j = 0; j < clusters[i].size(); j++) {
         entCluster.push_back(entList[clusters[i][j]]);
      }
      std::sort(entCluster.begin(), entCluster.end());
      for (uint32_t j = 0; j < entCluster.size(); j++) {
         printf("cl  %20s %4.2fA %4u %2u %s\n", entCluster[j].intname.c_str(),
               entCluster[j].resolution, entCluster[j].year,
               entCluster[j].numNonHohLig, entCluster[j].title.c_str());
      }
   }

   for (uint32_t i = 0; i < clusters.size(); i++) {
         const std::vector<float> dists = getMedoidDist(clusters[i], dm);
         for (uint32_t j = 0; j < clusters[i].size(); j++) {
            UqEntry &uq = entList[clusters[i][j]];
            /* add one to make node IDs associated with dendrogram */
            printf("%-15s %u %u %5.3f\n", uq.intname.c_str(), i, j + 1, dists[j]);
         }
      }

   delete dm;

   return 0;
}
