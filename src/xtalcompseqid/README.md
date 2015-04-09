For a list of PDB files it can compute three flat_file tables:
* chain connectivity - for all pairs of chains of the same PDB file, it lists the size of the interface (if > _threshold_)
* chain groups -  for all pairs of chains of the same PDB file, it indicates of they have (almost) identical sequences
* chain similarity - for all pairs of chains of different PDB files, it lists the sequence identity (if > _threshold_)

With these three tables one can generate a list of potential protein complexes with corresponding unbound structures (see xtalcompunbound) by using relational algrebra.

* Runs in parallel (using openmp) on all available cores
* Chunk size is optimized to avoid concurrent reads
