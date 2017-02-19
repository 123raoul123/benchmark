#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../test.h"
#include "zeropad_mul.h"
#include "fft/fiduccia.h"
#include "fft/split_radix_fft.h"
#include "fft/twisted_fft.h"
#include "fft/tangent_fft.h"

/******************************************************************
*
* TANGENT FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void zeropad_tangent(ring_t *r, const ring_t *x, const ring_t *y)
{ 
  double complex cplx_x[ZEROPAD];
  double complex cplx_y[ZEROPAD];
  double complex cplx_res[ZEROPAD];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
    cplx_x[i+REALDIM] = 0.0;
    cplx_y[i+REALDIM] = 0.0;
  }

  tangent_forward(cplx_x,ZEROPAD);
  tangent_forward(cplx_y,ZEROPAD);

  for (int i = 0; i < ZEROPAD; ++i)
  {
    cplx_res[i] = cplx_x[i] * cplx_y[i];
  }
  tangent_backward(cplx_res,ZEROPAD);

  for (int i = 0; i < REALDIM; ++i)
  {
    r->v[i] = (cplx_res[i] - cplx_res[i+REALDIM])/ZEROPAD;
  }
}

/******************************************************************
*
* SPLIT RADIX FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void zeropad_split_radix(ring_t *r, const ring_t *x, const ring_t *y)
{	
  double complex cplx_x[ZEROPAD];
  double complex cplx_y[ZEROPAD];
  double complex cplx_res[ZEROPAD];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
    cplx_x[i+REALDIM] = 0.0;
    cplx_y[i+REALDIM] = 0.0;
  }

  split_radix_recursive(cplx_x,ZEROPAD,0);
  split_radix_recursive(cplx_y,ZEROPAD,0);

  for (int i = 0; i < ZEROPAD; ++i)
  {
    cplx_res[i] = cplx_x[i] * cplx_y[i];
  }
  split_radix_recursive_inverse(cplx_res,ZEROPAD,0);

  for (int i = 0; i < REALDIM; ++i)
  {
    r->v[i] = (cplx_res[i] - cplx_res[i+REALDIM])/ZEROPAD;
  }
}

/******************************************************************
*
* Twisted FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void zeropad_twisted(ring_t *r, const ring_t *x, const ring_t *y)
{ 
  double complex cplx_x[ZEROPAD];
  double complex cplx_y[ZEROPAD];
  double complex cplx_res[ZEROPAD];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
    cplx_x[i+REALDIM] = 0.0;
    cplx_y[i+REALDIM] = 0.0;
  }

  twisted_recursive(cplx_x,ZEROPAD,0);
  twisted_recursive(cplx_y,ZEROPAD,0);

  for (int i = 0; i < ZEROPAD; ++i)
  {
    cplx_res[i] = cplx_x[i] * cplx_y[i];
  }
  twisted_recursive_inverse(cplx_res,ZEROPAD,0);

  for (int i = 0; i < REALDIM; ++i)
  {
    r->v[i] = (cplx_res[i] - cplx_res[i+REALDIM])/ZEROPAD;
  }
}

/******************************************************************
*
* SCHOOLBOOK NEGACYCLIC FFT
*
******************************************************************/
void zeropad_fiduccia(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[ZEROPAD];
  double complex cplx_y[ZEROPAD];
  double complex cplx_res[ZEROPAD];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
    cplx_x[i+REALDIM] = 0.0;
    cplx_y[i+REALDIM] = 0.0;
  }

  recursive_phi(cplx_x,ZEROPAD,0,1);
  recursive_phi(cplx_y,ZEROPAD,0,1);

  for (int i = 0; i < ZEROPAD; ++i)
  {
    cplx_res[i] = cplx_x[i] * cplx_y[i];
  }
  inverse_phi(cplx_res,ZEROPAD,0,1);

  for (int i = 0; i < REALDIM; ++i)
  {
    r->v[i] = (cplx_res[i] - cplx_res[i+REALDIM])/ZEROPAD;
  }
}

/******************************************************************
*
*	NAIVE SCHOOLBOOK MULTIPLICATION
*
******************************************************************/
/* Very simple schoolbook multiplication. Works. */
void schoolbook(ring_t *r, const ring_t *x, const ring_t *y)
{
  int i,j;
  for(i=0;i<REALDIM;i++)
    r->v[i] = 0;

  for(i=0;i<REALDIM;i++)
  {
    for(j=0;j<REALDIM;j++)
    {
      if(i+j < REALDIM)
        r->v[i+j] += x->v[i] * y->v[j];
      else
        r->v[i+j-REALDIM] -= x->v[i] * y->v[j];
    }
  }
}

