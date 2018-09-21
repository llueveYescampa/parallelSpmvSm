#include <stdlib.h>
#include <mpi.h>
#include <string.h>

void allocateSingleSharedVector(void **v, size_t size, const int *n, MPI_Win *sm_win, MPI_Comm *sm_comm) 
{

    int sharedRank;
    MPI_Comm_rank(*sm_comm,&sharedRank);

    MPI_Info myInfo;
    MPI_Info_create(&myInfo);
    MPI_Info_set(myInfo,"alloc_shared_noncontig", "false");
    
    if (sharedRank == 0) {
        MPI_Win_allocate_shared((MPI_Aint) (*n)*size,size,MPI_INFO_NULL,*sm_comm, v,sm_win);
    } else {
        MPI_Win_allocate_shared((MPI_Aint) 0,                size,MPI_INFO_NULL,*sm_comm, v,sm_win);
    } // end if //
        
    MPI_Aint sz;
    int dispUnit;

    MPI_Win_shared_query(*sm_win, MPI_PROC_NULL, &sz,&dispUnit,v);
    MPI_Win_lock_all(0,*sm_win);
    
    if (sharedRank == 0) {
        memset(*v, 0, (*n)*size  );
    } // end if //
    MPI_Win_sync(*sm_win);
    MPI_Barrier(*sm_comm);
} // end of allocateSingleSharedVector() //
