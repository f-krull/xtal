#ifndef CONFIGCOMPUNBOUND_H_
#define CONFIGCOMPUNBOUND_H_

#include <stdint.h>

class ConfigCompUnbound {
public:

   ConfigCompUnbound() {
      /* two chains are connected if their interface is above this threshold */
      minIntSize = 25;

      /* check if pair has similar sequence */
      chainMatchMinSeqId = 0.8;

      /* chainmatch ignores everything below */
      minSeqLen = 7; // according to Wikipedia MHC peptides have a size of 8-10

      /* all pairs of interfaces within a cluster have a distance < uqMaxDist */
      uqClMaxDist = 1 - 0.4f; /* 1-x; lower x -> less clusters */
      // cl max: 1 - 0.1 -> 2000 clusters


      /* no interface with (distance < u1MinIntDist) to B1:B2 is allowed within U1:^U1 */
      u1MinIntDist = 1 - 0.5f;

      /* how close to u1 are b2 chains allowed to be? */
      u1b2SeqIdMax = 0.7;

      /* allow U1 -> B2 clashes for oligomers */
      fixOligomerConflicts = true;

      /* of two chains have a seq ID above x they can be a homoDimer; heterodimer otherwise */
      u1HomodimerSeqidMin = 0.93;
      /* corner case:
       *  2cfh1 C A 2bjn1 A | ChainMatchExpansion-U1oligofix: B(u) A(u) 0.937063
       * */

      /* U1 -> B2 at most a fraction of x ca atoms are allowed to clash */
      u1ToB2MaxConflictRatio = 0.25;

      /* ca distances below x are considered as conflict */
      u1ToB2CaClashDist = 5.0;

      /* min distance of two residues in order be detected as part of interface (chainmatch edges only) */
      cmIntfDefinition = new IntfDefAllAtom(5.5);

      /* a matched edge is of if interface score of (b1,bi)(u1,ui) > x */
      cmMinInfSim = 0.3;

      /* DNA is not allowed to be closer than x to any interface residue */
      minDnaDist = 6.0f;

      /* cofactor are detected within interface if they have x neighboring interface residues within a cutoff of y A  */
      minCofactorDistB = 5.5f;
      minCofactorNeighborsB = 3;
      minCofactorDistU = 6.0f;
      minCofactorNeighborsU = 2;
      /* neighborhood of cofactor used for matching */
      cofNeighborDist = 10.0f;

      /* interface cutoff used for iRMSD and superimposition (compatibility to BM4) */
      superimoseIntfDefinition = new IntfDefMendez2003();

      /* how many residues need to be present in unbound */
      minB1U1alignedInterfaceRatio = 0.55f;

      maxU1IRmsd = 16.0f;

      gapMaxCtoNdist = 2.5f;

      cmxUCfixMaxResis = 30;
      cmxUCfixMaxS2bonds = 0;
      fixSmallPeptideClashes = false;
   }

   ~ConfigCompUnbound() {
      delete cmIntfDefinition;
      delete superimoseIntfDefinition;
   }


   uint32_t minIntSize;
   float chainMatchMinSeqId;
   uint32_t minSeqLen;
   float u1HomodimerSeqidMin;
   float uqClMaxDist;

   float u1b2SeqIdMax;
   float u1ToB2MaxConflictRatio;
   float u1ToB2CaClashDist;
   IntfDefAllAtom *cmIntfDefinition;
   float cmMinInfSim;
   float u1MinIntDist;
   float minDnaDist;
   float minCofactorDistB;
   uint32_t minCofactorNeighborsB;
   float minCofactorDistU;
   uint32_t minCofactorNeighborsU;
   float cofNeighborDist;

   IntfDefMendez2003 *superimoseIntfDefinition;
   float minB1U1alignedInterfaceRatio;
   float maxU1IRmsd;
   float gapMaxCtoNdist;

   bool fixOligomerConflicts;
   bool fixSmallPeptideClashes;
   uint32_t cmxUCfixMaxResis;
   uint32_t cmxUCfixMaxS2bonds;



private:
};



#endif /* CONFIGCOMPUNBOUND_H_ */
