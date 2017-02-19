#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "support.h"


#define omg_t(N,k) (cexp(2.0 * M_PI * I * (double)k / (double) N))


void twist(double complex *cplx_x,int n,int m,int lo)
{
  // printf("n = %d, m = %d, lo = %d\n",n,m,lo );
  int j = 1;
  for (int i = lo+1; i < lo+m; ++i)
  {
    // printf("i = %d, j = %d\n",i,j);
    cplx_x[i] = cplx_x[i] * omg_t(n,j);
    ++j;
  }
}

void untwist(double complex *cplx_x,int n,int m,int lo)
{
  // printf("n = %d, m = %d, lo = %d\n",n,m,lo );
  int j = 1;
  for (int i = lo+1; i < lo+m; ++i)
  {
    // printf("i = %d, j = %d\n",i,j);
    cplx_x[i] = cplx_x[i] * conj(omg_t(n,j));
    ++j;
  }
}