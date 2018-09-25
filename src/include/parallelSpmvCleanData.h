    if (numberOfNodes>1) {
        MPI_Win_free(&smWin_recvCount);
        MPI_Win_free(&smWin_sendCount);

        for (int node=0; node<numberOfNodes; ++node) {
            MPI_Win_sync(  smWin_compressedVec[node]  );
        } // end for //
        MPI_Barrier(sm_comm);
        
        for (int node=0; node<numberOfNodes; ++node) {
            MPI_Win_unlock_all( smWin_compressedVec[node] );
        } // end for //
        
        
        for (int node=0; node<numberOfNodes; ++node) {
            MPI_Win_free(  smWin_compressedVec+node  );
            MPI_Win_free(  smWin_sendColumns+node  );
        } // end for //
        free(smWin_compressedVec);
        free(smWin_sendColumns);        
        free(sendColumns);
        free(compressedVec);
    } // end if //
    

    free(requestS);
    free(requestR);

    free(row_ptr);
    free(col_idx);
    free(val);
    free(row_ptr_off);
    free(col_idx_off);
    free(val_off);
    
    MPI_Win_unlock_all(sm_win);
    MPI_Win_free(&sm_win);

    MPI_Win_unlock_all(smWin_v_off_nodal);
    MPI_Win_free(&smWin_v_off_nodal);
    
