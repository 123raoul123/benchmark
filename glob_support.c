#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "test.h"

/******************************************************************
*
* SUPPORT CODE
*
******************************************************************/
void print_complex(const double complex *a, int N)
{
    for(int i=0;i<N;++i)
      printf("cplxpoly[%d] = %f + i * %f\n",i,creal(a[i]),cimag(a[i]));
    printf("\n");
}

void print_double(const cplx *x,int N)
{
  for (int i = 0; i < N; ++i)
    printf("cplxpoly[%d] = %f + i * %f\n",i,x->real[i],x->imag[i]);
  printf("\n");
}

/******************************************************************
*
* CONVERSION
*
// ******************************************************************/
// void to_complex(const ring_t *x, const ring_t *y, cplx *cplx_x, cplx *cplx_y)
// {
//   int j = CPLXDIM;
//   for (int i = 0; i < CPLXDIM; ++i)
//   {
//     cplx_x->real[i] = x->v[i];
//     cplx_y->real[i] = y->v[i];
//     cplx_x->imag[i] = x->v[j];
//     cplx_y->imag[i] = y->v[j];
//     ++j;
//   }
// }

void to_real(const double complex *cplx_x, ring_t *x)
{
  for(int i=0;i<CPLXDIM;++i){
    x->v[i] = round(creal(cplx_x[i]));
    x->v[i+CPLXDIM] = round(cimag(cplx_x[i]));
  }
}

