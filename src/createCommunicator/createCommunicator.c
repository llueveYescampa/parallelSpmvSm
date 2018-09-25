#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "parallelSpmv.h"
#include "real.h"


int createColIdxMap(int **b, MPI_Win *sm_win, int *a, const int *n);

void createCommunicator(int *nColsOff,
                        int **recvCount,    MPI_Win *smWin_recvCount,
                        int **sendCount,    MPI_Win *smWin_sendCount,
                        int ***sendColumns, MPI_Win **smWin_sendColumns,
                        int *col_idx_off,
                        const int *off_proc_nnz,
                        const int *rowsPerNode,
                        real ***compressedVec, MPI_Win **smWin_compressedVec,
                        const int *nnodes,
                        int *countR,
                        int *countS ,
                        MPI_Request **requestR,
                        MPI_Request **requestS
                        )
{
    int  worldRank;
    MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);

    // creating an intranode communicator 
    MPI_Comm sm_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED, 0,MPI_INFO_NULL, &sm_comm);
    int sharedRank,sharedSize;
    MPI_Comm_rank(sm_comm,&sharedRank);
    MPI_Comm_size(sm_comm,&sharedSize);
    // creating an intranode communicator 

    // creating a communicator for the lead processes in each node
    MPI_Comm nodeComm;
    MPI_Comm_split(MPI_COMM_WORLD, sharedRank, (worldRank/sharedSize), &nodeComm );
    int nodeNumber; 
    MPI_Comm_rank(nodeComm,&nodeNumber);
    // creating a communicator for the lead processes in each node


    int *off_node_column_map=NULL;
    MPI_Win smWin_off_node_column_map;

    int off_node_nnz=0;
    MPI_Allreduce( off_proc_nnz, &off_node_nnz,1,MPI_INT,MPI_SUM,sm_comm);

    // creating the off_node_column_map as a shared array.
    if (off_node_nnz)  {
        // nColsOff is the number of off-node columns per node
        // only to be created if the node has off-node non-zeros
        *nColsOff=  createColIdxMap(&off_node_column_map,&smWin_off_node_column_map,col_idx_off,off_proc_nnz ); 
    } // end if //
    
    // creating the firstColumn array
    int *firstColumnArray=NULL;
    MPI_Win smWin_firstColumnArray;
    int size=(*nnodes+1);
    allocateSingleSharedVector((void **) &firstColumnArray,sizeof(int),  &size, &smWin_firstColumnArray, &sm_comm);
    //firstColumnArray  = (int *) malloc( (*nnodes+1)* sizeof(int)); 
    firstColumnArray[0]=0;
    if (sharedRank==0) {
        MPI_Allgather(rowsPerNode, 1,MPI_INT,&firstColumnArray[1], 1,MPI_INT,nodeComm);
        for (int i=1; i <= *nnodes; ++i) {
            firstColumnArray[i] += firstColumnArray[i-1];
        } // enf for //
        // end of creating the firstColumn array
    } // end if //
    MPI_Win_sync(smWin_firstColumnArray);
    MPI_Barrier(sm_comm);
    MPI_Win_unlock_all(smWin_firstColumnArray);
    // end of creating the firstColumn array


    // request arrays for communication  - two per node (send/recv)
    MPI_Request *myReqS=NULL,*myReqR=NULL;
    int **reciveColumns=NULL; 

    allocateSingleSharedVector( (void **) recvCount, sizeof(int),nnodes,  smWin_recvCount, &sm_comm) ;
    allocateSingleSharedVector( (void **) sendCount, sizeof(int),nnodes,  smWin_sendCount, &sm_comm) ;
    
    if (sharedRank==0) {
        //(*recvCount) = (int *) calloc( (unsigned int) *nnodes,sizeof(int)); 
        //(*sendCount) = (int *) calloc(*nnodes,sizeof(int));

        // finding rank holding each off_node column
        // and establishing the receive count arrays
        // process: is the rank
        // recvCount[process]: how many to receive from that rank 
        for (int i=0; i<*nColsOff; ++i) {
            int node=0;
            while( off_node_column_map[i] < firstColumnArray[node] ||  off_node_column_map[i] >= firstColumnArray[node+1]  ) {
                ++node;
            }// end while //
            ++(*recvCount)[node];
        } // end for //

        // establishing the send count arrays from the receive count arrays
        // node: is the node
        // sendCount[node]: how many to send to node node
        myReqS = (MPI_Request *) malloc(*nnodes*sizeof(MPI_Request));
        myReqR = (MPI_Request *) malloc(*nnodes*sizeof(MPI_Request));
            
        for (int node=0; node < *nnodes; ++node ) {
            MPI_Isend((*recvCount+node),1,MPI_INT,node, 123,nodeComm,&myReqS[node]);
            MPI_Irecv((*sendCount+node),1,MPI_INT,node, 123,nodeComm,&myReqR[node]);
        } // end for //

        // Crerating a 2d-array capable to hold rows of
        // independent size to store the lists of columns 
        // this rank need to receive.
        reciveColumns = (int **) malloc(*nnodes*sizeof(int *));
        for (int node=0; node<*nnodes; ++node){
            reciveColumns[node] = (int *) malloc((*recvCount)[node]*sizeof(int ));
        } // end for //
            
        // filling the reciveColumns arrays and 
        // adjustring it to make sendColumns relative to the current process
        for (int node=0, k=0; node<*nnodes; ++node){
            for (int i=0; i<(*recvCount)[node]; ++i, ++k){
                reciveColumns[node][i] = off_node_column_map[k]  - firstColumnArray[node];
            } // end for //
        } // end for //
        
        MPI_Waitall(*nnodes,myReqR,MPI_STATUS_IGNORE);
        MPI_Waitall(*nnodes,myReqS,MPI_STATUS_IGNORE);
    } // end if //


    // Crerating a 2d-array capable to hold rows of
    // independent size to store the lists of columns 
    // this rank need to send.
    (*sendColumns)   = (int **) malloc(*nnodes * sizeof(int *));
    (*smWin_sendColumns)   = (MPI_Win *) malloc(*nnodes * sizeof(MPI_Win *));
    
    for (int node=0; node<*nnodes; ++node) {
        allocateSingleSharedVector( (void **)  (*sendColumns)+node ,sizeof(int),  (*sendCount)+node ,  (*smWin_sendColumns)+node,  &sm_comm  ) ;
    } // end for //
    // end of establishing the send/receive count arrays

    if (sharedRank==0) {

        // Communicating the receive lists to sending processes to create the send lists
        for (int node=0; node <*nnodes; ++node) {
            if ((*recvCount)[node] > 0) {
                MPI_Isend(reciveColumns[node],(*recvCount)[node],MPI_INT,node, 321,nodeComm,&myReqS[node] );
            } // end if //   
            
            if ((*sendCount)[node] > 0) {
                MPI_Irecv(  (*sendColumns)[node],(*sendCount)[node],MPI_INT,node, 321,nodeComm,&myReqR[node] );
            } // end if //   
        } // end for //
        MPI_Waitall(*nnodes,myReqR,MPI_STATUS_IGNORE);
        MPI_Waitall(*nnodes,myReqS,MPI_STATUS_IGNORE);
        // end of Communicating the receive lists to sending processes to create the send lists
        free(myReqS);
        free(myReqR);

        for (int node=0; node<*nnodes; ++node){
            free(reciveColumns[node]);
        } // end for /
        free(reciveColumns);
    } // end if //
    MPI_Win_free(&smWin_firstColumnArray);

    // finaly, compressedVec is an 1D-array holding values on x_ptr to send to 
    // other process to build the x_off_ptr. The compressedVec could be as large
    // as the number of columns(rows) of this rank. It is used during the solution
    // to in the Isend/Irecv phase of the Spmv solution.
    (*compressedVec) = (real **) malloc(*nnodes*sizeof(real *));
    (*smWin_compressedVec)   = (MPI_Win *) malloc(*nnodes * sizeof(MPI_Win *));
    
    
    for (int node=0; node<*nnodes; ++node) {
        //(*compressedVec)[node] = (real *) malloc((*sendCount)[node]*sizeof(real ));
        allocateSingleSharedVector( (void **) (*compressedVec)+node ,sizeof(real),  (*sendCount)+node ,  (*smWin_compressedVec)+node, &sm_comm  ) ;
    } // end for //


    MPI_Win_sync(*smWin_recvCount);
    MPI_Win_sync(*smWin_sendCount);
    for (int node=0; node<*nnodes; ++node) {
        MPI_Win_sync(  (*smWin_sendColumns)[node]  );
    } // end for //
    
    MPI_Barrier(sm_comm);
    MPI_Win_unlock_all(*smWin_recvCount);
    MPI_Win_unlock_all(*smWin_sendCount);
    for (int node=0; node<*nnodes; ++node) {
        MPI_Win_unlock_all((*smWin_sendColumns)[node]);
    } // end for //
    
    if (off_node_nnz)  {
        MPI_Win_free(&smWin_off_node_column_map);
    } // end if //

//////////////////

    if (sharedRank==0 )   {
        for (int node=0; node<*nnodes; ++node) {
            if ((*recvCount)[node] > 0 ) ++*countR;
            if ((*sendCount)[node] > 0 ) ++*countS;
        } // end for //
        *requestS = (MPI_Request *) malloc( *countS*sizeof(MPI_Request));
        *requestR = (MPI_Request *) malloc( *countR*sizeof(MPI_Request));
    } // end if //


/////////////////////

} // end of createCommunicator() //




























/*
    if (sharedRank  % sharedSize ==  0) {
        printf("sharedRank: %d, ---> ", sharedRank);
        for (int node=0; node< *nnodes; ++node) {
            for (int i=0; i<(*sendCount)[node]; ++i) {
                printf("%4d", (*sendColumns)[node][i]  );
            }
            printf("\n");
        }
    }
    MPI_Finalize();
    exit(0);
*/
/*
    if (worldRank==2) {
        printf("WosdsdsrldRank: %d, ---> ", worldRank);
        for (int node=0; node< *nnodes; ++node) {
            for (int i=0; i<(*sendCount)[node]; ++i) {
                printf("%4d", (*sendColumns)[node][i]  );
            }
            printf("\n");
        }
    }
    MPI_Finalize();
    exit(0);

*/


/*
    printf("Rank: %d,  off_proc_nnz %d   --> ", worldRank, *off_proc_nnz ) ;
    for (int i=0; i< *off_proc_nnz; ++i){
        printf("%4d", col_idx_off[i]);
    }
    printf("\n");
    MPI_Finalize();
    exit(0);
*/   
/*
 /////////////////////////////////////////
    printf("RRRRrank: %d,  nColsOff %d --> ", worldRank,*nColsOff ) ;
    for (int i=0; i< *nColsOff; ++i){
        printf("%4d", off_node_column_map[i]);
    }
    printf("\n");
    MPI_Finalize();
    exit(0);
 /////////////////////////////////////////
*/

