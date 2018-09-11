#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "real.h"

#include "parallelSpmv.h"

#define REP 1000


int main(int argc, char *argv[]) 
{
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    int worldRank, worldSize;
    MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
    MPI_Comm_size(MPI_COMM_WORLD,&worldSize);

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
    int nodeNumber,numberOfNodes;
    MPI_Comm_rank(nodeComm,&nodeNumber);
    MPI_Comm_size(nodeComm,&numberOfNodes);
    // creating a communicator for the lead processes in each node

    const int root=0;

    
    if (worldRank == 0) {
        printf("num of nodes used: %d, worldSize: %d, sharedSize: %d\n", numberOfNodes,worldSize,sharedSize);
        //printf("worldRank %d, I am in node: %d out of %d nodes\n",worldRank,nodeNumber,numberOfNodes );
    } // end if //

    #include "parallelSpmvData.h"

    // verifing number of input parameters //
    char exists='t';
    char checkSol='f';
    if (worldRank == root) {
        if (argc < 3 ) {
            printf("Use: %s  Matrix_filename InputVector_filename  [SolutionVector_filename]  \n", argv[0]);     
            exists='f';
        } // endif //
        
        FILE *fh;
        // testing if matrix file exists
        if((fh = fopen(argv[1], "rb")  )   == NULL) {
            printf("No matrix file found.\n");
            exists='f';
        } // end if //
        
        // testing if input file exists
        if((fh = fopen(argv[2], "rb")  )   == NULL) {
            printf("No input vector file found.\n");
            exists='f';
        } // end if //

        // testing if output file exists
        if (argc  >3 ) {
            if((fh = fopen(argv[3], "rb")  )   == NULL) {
                printf("No output vector file found.\n");
                exists='f';
            } else {
                checkSol='t';
            } // end if //
        } // end if //
    } // end if //
    MPI_Bcast(&exists,  1,MPI_CHAR,root,MPI_COMM_WORLD);
    if (exists == 'f') {
       if (worldRank == root) printf("Quitting.....\n");
        MPI_Finalize();
        exit(-1);
    } // end if //
    MPI_Bcast(&checkSol,1,MPI_CHAR,root,MPI_COMM_WORLD);

    int n_global,nnz_global;
    int rowsPerProc,rowsPerNode;
    int off_proc_nnz=0;
    
    reader(&n_global,&nnz_global, &rowsPerProc, 
           &off_proc_nnz,
           &row_ptr,&col_idx,&val,
           &row_ptr_off,&col_idx_off,&val_off,
           argv[1], root);
           
    // determining the number of rows per node//           
    MPI_Allreduce(&rowsPerProc,&rowsPerNode,1,MPI_INT,MPI_SUM,sm_comm );
    // end of determining the number of rows per node//           

    // nColsOff is the number of off-node columns per node
    int nColsOff=0;
    if (worldSize>1) {
        createCommunicator(&nColsOff, &recvCount,&smWin_recvCount, &sendCount,&smWin_sendCount, &sendColumns,&smWin_sendColumns, col_idx_off, &off_proc_nnz, &rowsPerNode,&compressedVec, &smWin_compressedVec, &numberOfNodes);
    } // end if //
        
    // ready to start //    

///////////////////////////////////////////////////
////  shared memory

    MPI_Win sm_win;
    real *v=NULL; 
    real *v_nodal=NULL; 

    MPI_Win smWin_v_off_nodal;
    real *v_off_nodal=NULL; 
    
    
    allocateSharedVector(&v,  &v_nodal,     &rowsPerProc, &sm_win,     &sm_comm );
    allocateSingleSharedVector((void **) &v_off_nodal,sizeof(real), &nColsOff,    &smWin_v_off_nodal, &sm_comm );
    
    //printf("worldRank: %d, nColsOff: %d\n", worldRank,nColsOff);
    
    
    //v_off_nodal = (real *) malloc((nColsOff)*sizeof(real));

    // reading input vector
    vectorReader(v, &rowsPerProc, argv[2]);
    MPI_Win_sync(sm_win);
    MPI_Barrier(sm_comm);  
    
///////////////////////////////////////////////////

    
    real *w=NULL;
    w     = (real *) malloc(rowsPerProc*sizeof(real)); 


    int countR=0, countS=0;
    if (sharedRank==0 && worldSize>1)   {
        for (int node=0; node<numberOfNodes; ++node) {
            if (recvCount[node] > 0 ) ++countR;
            if (sendCount[node] > 0 ) ++countS;
        } // end for //
        requestS = (MPI_Request *) malloc( countS*sizeof(MPI_Request));
        requestR = (MPI_Request *) malloc( countR*sizeof(MPI_Request));
    } // end if //
    

    // Timing should begin here//
    double elapsed_time;
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();
    
    for (int t=0; t<REP; ++t) {
        // cleaning solution vector //
        for(int i=0; i<rowsPerProc; ++i) w[i] = 0.0;

        
            if (worldSize>1) {
                startComunication(v_nodal,v_off_nodal,compressedVec,recvCount, sendCount, sendColumns, requestS,requestR,&numberOfNodes, &nodeComm, &sharedRank);
            }  // end if // 

        // solving the on_proc part while comunication is taken place.
        spmv(w,val,v_nodal, row_ptr,col_idx,rowsPerProc);

        
        // waitting for the comunication to finish
        if (sharedRank == 0 && worldSize>1 ) {
            MPI_Waitall(countR, requestR,MPI_STATUS_IGNORE);
            MPI_Waitall(countS, requestS,MPI_STATUS_IGNORE);
        } // end if //
        

        MPI_Win_sync(smWin_v_off_nodal);
        MPI_Barrier(sm_comm);
        // now is time to solve the off_proc part
        if (off_proc_nnz > 0) {
            spmv(w,val_off,v_off_nodal, row_ptr_off,col_idx_off,rowsPerProc);
        } // end if//

        MPI_Barrier(MPI_COMM_WORLD);
    } // end for //
    // Timing should end here//
    elapsed_time += MPI_Wtime();



    if (worldRank == root) {
        printf("---> Time taken by %d processes: %g seconds\n",worldSize, elapsed_time);
    } // end if //

    
    if (checkSol=='t') {
        real *sol=NULL;
        sol     = (real *) malloc((rowsPerProc)*sizeof(real)); 
        // reading input vector
        vectorReader(sol, &rowsPerProc, argv[3]);
        
        int row=0;
        const real tolerance=1.0e-08;
        real error;
        do {
            error =  fabs(sol[row] - w[row]); //  /fabs(sol[row]);
            if ( error > tolerance ) break;
            ++row;
        } while (row < rowsPerProc); // end do-while //
        
        if (row == rowsPerProc) {
            printf("Solution match in rank %d\n",worldRank);
        } else {    
            printf("For Matrix %s, solution does not match at element %d in rank %d   %20.13e   -->  %20.13e  error -> %20.13e, tolerance: %20.13e \n", 
            argv[1], (row+1),worldRank, sol[row], w[row], error , tolerance  );
        } // end if //
        free(sol);    
    } // end if //
    
 

    free(w);
    
    #include "parallelSpmvCleanData.h" 
    MPI_Finalize();
    return 0;    
} // end main() //


/*



    if (sharedRank==0)   {
        printf("wordrank: %d -->", worldRank); 
        for (int node=0; node < numberOfNodes;++node ) {
            for (int i=0;  i<sendCount[node]; ++i) {
                printf("%6.2f",compressedVec[node][i]  );
            } // end for //
            printf("\n");
        } // end for //
    } // end if //
    MPI_Finalize(); exit(0);



    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    if (sharedRank==0) {
        printf("%s, rank: %d, local rank: %d, -->  !|", processor_name,worldRank,sharedRank ) ;
        for (int node=0; node<numberOfNodes; ++node) {
            for (int j=0; j<sendCount[node]; ++j) {
                printf("%4d", sendColumns[node][j] );
            } //
            printf("|\t\t");
        } //
        printf("\n");
    } // end if


    MPI_Finalize();
    exit(0);

*/




/*    

    printf("rank: (%d) ---> ", worldRank);
    for (int j=0; j<n; ++j) {
        printf("%f  ", w[j] );
    } //
    printf("\n");



    printf("rank: (%d) ---> ", worldRank);
    for (int process=0; process<worldSize; ++process) {
        for (int j=0; j<sendCount[process]; ++j) {
            printf("%4d", sendColumns[process][j] );
        } //
        printf("|\t\t");
    } //
    printf("\n");
    MPI_Finalize();
    exit(0);




    /////////////////////////////////////////////////////////////       

    int node_nnz=0;
    for (int row=0; row<n; ++row) {
        node_nnz += row_ptr[row+1] - row_ptr[row];
    }
    
    printf("rank: %d, local rank: %d, on_node_nnz %d, -->", worldRank,sharedRank,node_nnz);
    for (int i=0; i<node_nnz; ++i) {
        printf("%4.0f",val[i] );
    } // end for //


    printf(",  off_node_nnz %d -->", off_node_nnz );
    for (int i=0; i<off_node_nnz; ++i) {
        printf("%4.0f",val_off[i] );
    } // end for //
    printf("\n");
    
    MPI_Finalize();
    exit(0);
    /////////////////////////////////////////////////////////////       
               


    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    printf("%s, rank: %d, local rank: %d, -->  !|", processor_name,worldRank,sharedRank ) ;
    for (int process=0; process<worldSize; ++process) {
        for (int j=0; j<sendCount[process]; ++j) {
            printf("%4d", sendColumns[process][j] );
        } //
        printf("|\t\t");
    } //
    printf("\n");


    MPI_Finalize();
    exit(0);


    // testing the values of v
    int nodeRows;
    MPI_Allreduce(&n,&nodeRows, 1, MPI_INT, MPI_SUM,sm_comm);
    
    
    printf("rank: %d, sharedRank: %d --> ", worldRank, sharedRank);
    
    
    for (int row=0; row<nodeRows; ++row) {
        printf("%12.5f,",  v_nodal[row]);
    }
    printf("\n");
    MPI_Finalize();
    exit(0);
    // testing the values of v

*/    

