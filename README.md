# parallelSpmvSm
This corresponds to the shared memory "MPI+MPI" version

Version 1
    First shared memory version. Each process potentially communicates with every other process in-node and out-node.

Version 2
    In version 2 Only one process per node (Rank 0) communicates with one process (Rank 0) in other nodes sending all the data that node needs.
    
Version 3
    The spmv version was modified to allow vectorization in inner loop. For the pgi compiler a new parameter was used to stop the check for dependency in the loop.
    The results for this version are good for Intel, mixed for pgi and no-good for gcc 7.3

