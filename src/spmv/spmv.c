#include "real.h"

void spmv(real *__restrict__  b, real *__restrict__ val, real *x, int *row_ptr, int *col_idx, int nRows, const real alpha, const real beta)
{
    int start, end;
    for (int i=0; i<nRows; ++i) {
        start = row_ptr[i];
        end = row_ptr[i+1];
        real dot = (real) 0.0;
        for (int j=start; j<end; ++j) {
            dot += val[j] * x[col_idx[j]];
        } // end for //
        b[i] = beta*b[i] + alpha*dot;
    } // end for //
} // end of spmv() //

