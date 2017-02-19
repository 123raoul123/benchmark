#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../test.h"
#include "../glob_support.h"
#include "lut_mul.h"
#include "fft/lut_negacyclic.h"
#include "fft/lut_split_radix.h"
#include "fft/lut_tangent.h"

void lut_init(){
  //PRECOMP TABLES FOR PRECOMP FFT
  init_table();
  //INIT LOOKUPTABLES NEGACYCLIC FFT
  init_negacyc();
  //INIT LOOKUPTABLES TANGENT FFT
  init_tangent();
}

/******************************************************************
*
* TANGENT FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void lut_tangent(ring_t *r, const ring_t *x, const ring_t *y)
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

  lut_tangent_forward(cplx_x);
  // print_complex(cplx_x,CPLXDIM);
  lut_tangent_forward(cplx_y);

  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }
  lut_tangent_backward(cplx_res);
  // print_complex(cplx_res,CPLXDIM);
  to_real(cplx_res,r);
}

/******************************************************************
*
* SPLIT RADIX PRECOMPUTED FFT MULTIPLICATION
*
******************************************************************/
void lut_split_radix(ring_t *r, const ring_t *x, const ring_t *y){
  cplx cplx_x,cplx_y,cplx_res;

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_x.real[i] = x->v[i];
    cplx_y.real[i] = y->v[i];
    cplx_x.imag[i] = x->v[j];
    cplx_y.imag[i] = y->v[j];
    ++j;
  }
  
  fft_precompsr_forward(&cplx_x);
  fft_precompsr_forward(&cplx_y);

  double a,b,c,d;
  for (int i = 0; i < CPLXDIM; ++i)
  {
      a = cplx_x.real[i];
      b = cplx_x.imag[i]; 
      c = cplx_y.real[i];
      d = cplx_y.imag[i];
      //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
      cplx_res.real[i] = ((a*c) - (b*d))/CPLXDIM;
      cplx_res.imag[i] = ((a*d) + (b*c))/CPLXDIM;
  }
  fft_precompsr_backward(&cplx_res,r);
}

/******************************************************************
*
* NEGACYCLIC FFT LOOK UP TABLE
*
******************************************************************/
void lut_negacyclic(ring_t *r, const ring_t *x, const ring_t *y)
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
  // printf("*****************START LUT FFT********************\n");
  recursive_phi_lut(cplx_x,CPLXDIM,0,0);
  // print_complex(cplx_x,CPLXDIM);
  recursive_phi_lut(cplx_y,CPLXDIM,0,0);


  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }

  inverse_phi_lut(cplx_res,CPLXDIM,0,0);

  to_real(cplx_res,r);
}