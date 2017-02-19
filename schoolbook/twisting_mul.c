#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../test.h"
#include "../glob_support.h"
#include "support.h"
#include "twisting_mul.h"
#include "fft/fiduccia.h"
#include "fft/split_radix_fft.h"
#include "fft/twisted_fft.h"
#include "fft/tangent_fft.h"

/******************************************************************
*
* TANGENT FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void twist_tangent(ring_t *r, const ring_t *x, const ring_t *y)
{
// { printf("*********************STARTING TANGENT***********************\n");
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_x[i] = x->v[i] +I*x->v[j];
    cplx_y[i] = y->v[i] +I*y->v[j];
    ++j;
  }

  twist(cplx_x,ZEROPAD,CPLXDIM,0);
  twist(cplx_y,ZEROPAD,CPLXDIM,0);

  tangent_forward(cplx_x,CPLXDIM);
  // print_complex(cplx_x,CPLXDIM);
  tangent_forward(cplx_y,CPLXDIM);

  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }
  tangent_backward(cplx_res,CPLXDIM);
  untwist(cplx_res,ZEROPAD,CPLXDIM,0);
  // print_complex(cplx_res,CPLXDIM);
  to_real(cplx_res,r);
}

/******************************************************************
*
* SPLIT RADIX FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void twist_split_radix(ring_t *r, const ring_t *x, const ring_t *y)
{	
  // printf("*********************STARTING SPLIT RADIX***********************\n");
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_x[i] = x->v[i] +I*x->v[j];
    cplx_y[i] = y->v[i] +I*y->v[j];
    ++j;
  }
  twist(cplx_x,ZEROPAD,CPLXDIM,0);
  twist(cplx_y,ZEROPAD,CPLXDIM,0);

  split_radix_recursive(cplx_x,CPLXDIM,0);
  split_radix_recursive(cplx_y,CPLXDIM,0);

  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }
  split_radix_recursive_inverse(cplx_res,CPLXDIM,0);
  untwist(cplx_res,ZEROPAD,CPLXDIM,0);
  // print_complex(cplx_res,CPLXDIM);
  to_real(cplx_res,r);
}

/******************************************************************
*
* Twisted FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void twist_twisted(ring_t *r, const ring_t *x, const ring_t *y)
{ 
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_x[i] = x->v[i] +I*x->v[j];
    cplx_y[i] = y->v[i] +I*y->v[j];
    ++j;
  }
  twist(cplx_x,ZEROPAD,CPLXDIM,0);
  twist(cplx_y,ZEROPAD,CPLXDIM,0);

  twisted_recursive(cplx_x,CPLXDIM,0);
  twisted_recursive(cplx_y,CPLXDIM,0);

  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }
  twisted_recursive_inverse(cplx_res,CPLXDIM,0);
  untwist(cplx_res,ZEROPAD,CPLXDIM,0);
  to_real(cplx_res,r);
}

/******************************************************************
*
* FIDUCCIA FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void twist_fiduccia(ring_t *r, const ring_t *x, const ring_t *y)
{ 
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_x[i] = x->v[i] +I*x->v[j];
    cplx_y[i] = y->v[i] +I*y->v[j];
    ++j;
  }
  twist(cplx_x,ZEROPAD,CPLXDIM,0);
  twist(cplx_y,ZEROPAD,CPLXDIM,0);

  recursive_phi(cplx_x,CPLXDIM,0,1);
  recursive_phi(cplx_y,CPLXDIM,0,1);

  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }
  inverse_phi(cplx_res,CPLXDIM,0,1);
  untwist(cplx_res,ZEROPAD,CPLXDIM,0);
  to_real(cplx_res,r);
}

/******************************************************************
*
* NEGACYCLIC FFT
*
******************************************************************/
void twist_negacyclic(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_x[i] = x->v[i] +I*x->v[j];
    cplx_y[i] = y->v[i] +I*y->v[j];
    ++j;
  }

  double complex root = I;
  root = csqrt(root);
  // print_complex(cplx_x,CPLXDIM);
  recursive_phi(cplx_x,CPLXDIM,0,root);
  // print_complex(cplx_x,CPLXDIM);
  recursive_phi(cplx_y,CPLXDIM,0,root);


  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }

  inverse_phi(cplx_res,CPLXDIM,0,root);

  to_real(cplx_res,r);
}
