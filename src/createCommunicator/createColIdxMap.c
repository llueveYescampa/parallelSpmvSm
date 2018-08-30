// this function permits eliminate repeated element of an integer array
// it also creates a new array with the correct size to hold only unique element
// the function retuens the size of the new array.

#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "parallelSpmv.h"

int cmpfunc (const void *a, const void *b) {
 // this function is to be used with the library qsort()
   return ( *(int*)a - *(int*)b );
} // end of cmpfunc //


int createColIdxMap(int **b, MPI_Win *sm_win, int *a, const int *n) 
{

    MPI_Comm sm_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED, 0,MPI_INFO_NULL, &sm_comm);

    int sharedRank,sharedSize;
    MPI_Comm_rank(sm_comm,&sharedRank);
    MPI_Comm_size(sm_comm,&sharedSize);

    int *temp=NULL;
    int *temp2=NULL;
    int count=0;
    if (*n > 0) {
        
        temp = (int *)  malloc( *n * sizeof(int));
        memcpy(temp, a, *n*sizeof(int));
        qsort(temp, *n, sizeof(int), cmpfunc);
        
        count=1;
        for (int i=1; i<*n;++i) {
            if (temp[i-1] != temp[i] ) {
                ++count;
            } // end if //
        } // end for //
        
        temp2 = (int *)  malloc( count * sizeof(int));
        
        temp2[0] = temp[0];
        for (int i=1,j=1; i<*n;++i) {
            if (temp[i-1] != temp[i] ) {
               temp2[j] = temp[i];
               ++j; 
            } // end if //
        } // end for //
        free(temp);
    } // end if //    
    
    int countNode=0;
    MPI_Reduce(&count,&countNode,1,MPI_INT,MPI_SUM,0,sm_comm);
    
    
    int *disp=NULL;
    int *recvCnt=NULL;
    if (sharedRank == 0) {
        disp    = (int *) calloc( (sharedSize+1) , sizeof(int));
        recvCnt = (int *) calloc( sharedSize , sizeof(int));
        temp    = (int *) calloc( countNode, sizeof(int));
    } // end if //
    
    MPI_Gather(&count,1,MPI_INT,&disp[1],1,MPI_INT,0,sm_comm );    
    MPI_Gather(&count,1,MPI_INT,recvCnt, 1,MPI_INT,0,sm_comm );    
    
    // creating disp //
    if (sharedRank == 0) {
        for (int i=0; i<sharedSize ; ++i) {
            disp[i+1] += disp[i];
        } // end for //
    } // end if //
    
    MPI_Gatherv( temp2,count,MPI_INT,temp,recvCnt,disp,MPI_INT,0,sm_comm);
    
    free(disp);
    free(recvCnt);
    free(temp2);
    

    if (sharedRank==0) {
        qsort(temp, countNode, sizeof(int), cmpfunc);
        count=1;
        for (int i=1; i<countNode  ;++i) {
            if (temp[i-1] != temp[i] ) {
                ++count;
            } // end if //
        } // end for //
    } // end if //
    MPI_Bcast(&count,1,MPI_INT,0,sm_comm);

    allocateSingleSharedVector( (void **)  b, sizeof(int),&count,  sm_win,  &sm_comm) ;
    
    if (sharedRank == 0) {
        (*b)[0] = temp[0];
        for (int i=1,j=1; i<  countNode;  ++i) {
            if (temp[i-1] != temp[i] ) {
               (*b)[j] = temp[i];
               ++j; 
            } // end if //
        } // end for //
    } // end if //
    MPI_Win_sync(*sm_win);
    MPI_Barrier(sm_comm);
    MPI_Win_unlock_all(*sm_win);
    if (sharedRank == 0)  free(temp);

    // compacting the sequence of input vector
    for (int j=0; j<*n; ++j) {
        // running a binary search in the b array
        int l = 0;
        int r = count-1;
        while ( l <= r ) {
            int m = (l+r)/2;
            if (  (*b)[m] < a[j] ) {
                l = m+1;
            } else if (  (*b)[m] > a[j] ) {
                r = m-1;
            } else {
                a[j] = m;
                break;
            } // end if //
        } // end while //
    }  // end for //
    return count;    
} // end of createColIdxMap() //
