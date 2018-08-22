
    // data for the on_proc solution
    int *row_ptr=NULL;
    int *col_idx=NULL;
    real *val=NULL;
    // end of data for the on_proc solution
    
    // data for the off_proc solution
    int *row_ptr_off=NULL;
    int *col_idx_off=NULL;
    real *val_off=NULL;
    // end of data for the off_proc solution

    // creatinng communicator data//
    int *recvCount=NULL; MPI_Win smWin_recvCount;
    int *sendCount=NULL;  MPI_Win smWin_sendCount;
    int **sendColumns=NULL; MPI_Win *smWin_sendColumns;
    // compressedVec is an 1D-array holding values on x_ptr to send to other process
    // to build the x_off_ptr. The compressedVec could be as
    // large as the number of columns(rows) of this rank
    real **compressedVec=NULL; MPI_Win *smWin_compressedVec;
    // end of creatinng communicator data//

    // request arrays for communication  - two per process (send/recv)
    MPI_Request *requestS=NULL, *requestR=NULL;
