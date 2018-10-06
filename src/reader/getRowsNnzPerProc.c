#include <math.h>
#include <mpi.h>

void getRowsNnzPerProc(int *rowsPP, int *nnzPP, const int *global_n, const int *global_nnz,  const int *row_Ptr) 
{
    int worldRank, worldSize;
    MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
    MPI_Comm_size(MPI_COMM_WORLD,&worldSize);

    float nnzIncre = (float ) *global_nnz/ (float) worldSize;
    float lookingFor=nnzIncre;
    int startRow=0, endRow;
    int partition=0;

    for (int row=0; row<*global_n; ++row) {    
        if ( (float) row_Ptr[row+1] >=  lookingFor ) { 
            // search for smallest difference
            if (fabs ( lookingFor - row_Ptr[row+1])  <= fabs ( lookingFor - row_Ptr[row])   ) {
                endRow = row;
            } else {
                endRow = row-1;
            } // end if //
            
            rowsPP[partition] = endRow-startRow+1;
            nnzPP[partition]  = row_Ptr[endRow+1] - row_Ptr[startRow];
             
            startRow = endRow+1;
            ++partition;
            if (partition < worldSize-1) {
               lookingFor += nnzIncre;
            } else {
                lookingFor=*global_nnz;
            } // end if //   
        } // end if // 
    } // end for //
} // end of getRowsPerProc //
/*
rowsPP[0] = 799356;
rowsPP[1] = 798021;
rowsPP[2] = 795534;
rowsPP[3] = 795654;
rowsPP[4] = 797769;
rowsPP[5] = 798276;
rowsPP[6] = 799338;
rowsPP[7] = 799338;
nnzPP[0] = 7930445;
nnzPP[1] = 6921748;
nnzPP[2] = 13426076;
nnzPP[3] = 9829236;
nnzPP[4] = 13269564;
nnzPP[5] = 9632820;
nnzPP[6] = 9921636;
nnzPP[7] = 8233430;
*/

