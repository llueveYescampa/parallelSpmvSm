#include <mpi.h>
#include "real.h"

void allocateSharedVector(real **v, real **v_nodal, const int *n, MPI_Win *sm_win, MPI_Comm *sm_comm ) 
{




    MPI_Info myInfo;
    MPI_Info_create(&myInfo);
    MPI_Info_set(myInfo,"alloc_shared_noncontig", "false");
    
    
    MPI_Win_allocate_shared((MPI_Aint) (*n)*sizeof(real),sizeof(real),myInfo,*sm_comm, v,sm_win);




    MPI_Aint sz;
    int dispUnit;
    
    MPI_Win_shared_query(*sm_win, MPI_PROC_NULL, &sz,&dispUnit,v_nodal);
    MPI_Win_lock_all(0,*sm_win);
    
    for (int row=0; row<*n; ++row ) {
        (*v)[row] = (real) 0.0;
    } // end for //    
    MPI_Win_sync(*sm_win);
    MPI_Barrier(*sm_comm);
    
} // end of allocateSharedVector() //    

