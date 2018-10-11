#include "real.h"
#define CEILING(i,j)       (((i)+(j)-1)/(j))

void reader(int *gn, int *gnnz, int *n,  int *off_proc_nnz, 
            int **rPtr,int **cIdx,real **v,int **rPtrO,int **cIdxO,real **vO,
            const char *matrixFile, const int root);

void vectorReader(real *v, const int *n, const char *vectorFile);

                        
void createCommunicator( int *nColsOff,
                        int **recvCount,    MPI_Win *smWin_recvCount,
                        int **sendCount,    MPI_Win *smWin_sendCount,
                        int ***sendColumns, MPI_Win **smWin_sendColumns,
                        int *col_idx_off,
                        const int *off_node_nnz,
                        const int *rowsPerNode,
                        real ***compressedVec, MPI_Win **smWin_compressedVec,
                        const int *nnodes,
                        int *countR,
                        int *countS,
                        MPI_Request **reqR,
                        MPI_Request **reqS,
                        int **ranks2Send,
                        int **ranks2Recv                         
                         );                        
                         

void startComunication(real *v,
                       real *v_off, 
                       real **compressedVec,
                       const int *recvCount,
                       const int *sendCount, 
                       int **sendColumns, 
                       MPI_Request *requestS,
                       MPI_Request *requestR,
                       const int *nnodes,
                       const int *sharedRank,
                       const int *sharedSize,
                       const int *ranks2Send,
                       const int *ranks2Recv
                        );

void spmv(real *b, real *__restrict__ val, real *x, int *row_ptr, int *col_idx, int nRows);                         
void allocateSharedVector(real **v, real **v_nodal,const int *n, MPI_Win *sm_win, MPI_Comm *sm_comm );
void allocateSingleSharedVector(void **v, size_t size, const int *n, MPI_Win *sm_win, MPI_Comm *sm_comm);
