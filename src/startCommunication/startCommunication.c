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
                       const int *sharedRank,
                       const int *sharedSize,
                       const int *ranks2Send,
                       const int *ranks2Recv
                        )
{
    for (int node=0,indx=0, s=0, r=0,tempR=0, tempS=0; node < *nnodes; indx += recvCount[node++] ) {
    
        // forming x_off_ptr array with contributions from each rank
        if (recvCount[node] > 0) {
            if ( *sharedRank == (*sharedSize - 1 - (tempR % *sharedSize)) ) {
                //printf("Node: %d, worldRank: %d, rank receiving: %d\n", node, worldRank, ranks2Recv[r]);
                MPI_Irecv(&x_off_ptr[indx],recvCount[node],MPI_MY_REAL,ranks2Recv[r], 231,MPI_COMM_WORLD,&requestR[r] );
                ++r;
            } // end if //
            ++tempR;
        } // end if //    

        // need to compress data to send inside compressedVec
        if (sendCount[node] > 0) {                   // need to send to process j
            for (int i=0;  i<sendCount[node]; ++i) {
                    compressedVec[node][i] = x_ptr[sendColumns[node][i]];
            } // end for //
            if ( *sharedRank == (tempS % *sharedSize) ) {
                //printf("Node: %d, worldRank: %d, rank sending: %d\n", node, worldRank,ranks2Send[s]);
                MPI_Isend(compressedVec[node],sendCount[node],MPI_MY_REAL, ranks2Send[s] , 231,MPI_COMM_WORLD,&requestS[s] );
                ++s;
            } // end if //
            ++tempS;                
        } // end if //
    } // end for //
} // end of startComunication() //
