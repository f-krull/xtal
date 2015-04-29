# xtaluniquecomp 

Computes the similarity of two protein docking complexes (*p1* and *p2*), where  p1 consists of the docking units *n1*, *m1* and *p2* consists of the two docking units *n2*, *m2*. The algorithm first computes sequence alignments between the docking units (*n1*->*n2* and *m1*->*m2*) and then analyzes the coverage of the residues within the interfaces. The cross-combination (*n1*->*m2* and *m1*->*n2*) is also considered.  
Instead of a pair of protein complexes a list of protein docking complexes can be provided. In this case, the algorithm computes all pairwise similarities and applies single-linkage clustering to all entries. 

## Usage
```
xtaluniquecomp p1 n1 m1 p2 n1 m2
```
* p1 path to PDB structure / protein docking complex with interface n1:m1
* n1 chains that belog to docking unit n1
* m1 chains that belog to docking unit m1
* p2 path to PDB structure / protein docking complex with interface n2:m2
* n2 chains that belog to docking unit n2
* m2 chains that belog to docking unit m2

## Example
```
# get test data
PDBURL=http://ftp.wwpdb.org//pub/pdb/data/structures/all/pdb/
wget -O- ${PDBURL}/pdb1eaw.ent.gz | gunzip > /tmp/1eaw1.pdb
wget -O- ${PDBURL}/pdb1fy8.ent.gz | gunzip > /tmp/1fy81.pdb

# compute similarity
src/xtaluniquecomp/xtaluniquecomp /tmp 1eaw1 A B 1fy81 E I

# the same similarity is computed by:
src/xtaluniquecomp/xtaluniquecomp /tmp 1eaw1 A B 1fy81 I E
src/xtaluniquecomp/xtaluniquecomp /tmp 1fy81 I E 1eaw1 A B 
```
