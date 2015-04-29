# xtalcompunbound

Takes the structure of two proteins a arguments. One structure (*B*) serves as the potential protein docking complex. The other structure (*U*) serves as the potential unbound structure. The algorithm verifies if the potential protein docking complex can be divided into two binding partners (*B1* and *B2*) such that the potential unbound structure *U* (partially) corresponds to  *B1*. 
To avoid multiple solutions, the algorithm additionally gets one initial chain from each binding partner and one initial chain from the unbound structure. 

## Usage
  ```
 xtalcompunbound pB cB1 cB2 pU cU1
  ```

* pB path to the structure of a potential protein docking complex 
* cB1 one chain of binding partner B1
* cB2 one chain of binding partner B2
* pU path to the structure of a potential unbound protein
* cU1 one chain corresponding to cB1

Requirements to find a solution:
* pB is a PDB structure with at leat two chains (cB1, and cB2)
* cB1 and cB2 are in contact
* pU is a PDB structure with at leat one chains (cU1)
* pB1 has high sequence identity to cU1
* pU has no chain with sequence identity to cB2

## Example
```
PDBURL=http://ftp.wwpdb.org//pub/pdb/data/structures/all/pdb/

wget -O- ${PDBURL}/pdb1akj.ent.gz | gunzip > /tmp/1akj.pdb
wget -O- ${PDBURL}/pdb2clr.ent.gz | gunzip > /tmp/2clr.pdb

# define similar cofactors 
echo "MN MG CO CA FE NI SR ZN CD HG" > /tmp/cof_groups.txt
# define cofactor not being considered
echo "HOH" > /tmp/cof_ignorelist.txt

src/xtalcompunbound/xtalcompunbound -debug=true /tmp/1akj.pdb A D /tmp/2clr.pdb D
```
