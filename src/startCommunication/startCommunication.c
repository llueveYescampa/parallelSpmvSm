#include <stdio.h>
#include <mpi.h>
#include "real.h"

void startComunication(real *x_ptr,
                       real *x_off_ptr,
                       real **compressedVec,
                       const int *recvCount,
                       const int *sendCount,  
                       int **sendColumns,
                       MPI_Request *requestS,
                       MPI_Request *requestR,
                       const int *nnodes,
                       MPI_Comm *nodeComm,
                       const int *sharedRank
                        )
{

    if (*sharedRank == 0) {
        for (int node=0,indx=0, s=0, r=0; node < *nnodes; indx += recvCount[node++] ) {
            // forming x_off_ptr array with contributions from each rank
            if (recvCount[node] > 0) {
                MPI_Irecv(&x_off_ptr[indx],recvCount[node],MPI_MY_REAL,node, 231,*nodeComm,&requestR[r++] );
            } // end if //    

            // need to compress data to send inside compressedVec
            if (sendCount[node] > 0) {                   // need to send to process j
                for (int i=0;  i<sendCount[node]; ++i) {
                        compressedVec[node][i] = x_ptr[sendColumns[node][i]];
                } // end for //
                MPI_Isend(compressedVec[node],sendCount[node],MPI_MY_REAL,node, 231,*nodeComm,&requestS[s++] );
            } // end if //
        } // end for //
    } // end if //
} // end of startComunication() //
