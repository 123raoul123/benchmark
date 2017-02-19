#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <immintrin.h>
#include "../test.h"
#include "vec_mul.h"
#include "fft/sr_vector.h"
#include "fft/fftw.h"
#include "fft/sr_vec_nonrec.h"
#include "fft/2_layer_negacyclic.h"
#include "fft/simd_negacyclic.h"
#include "fft/3_layer_negacyclic.h"
#include "fft/simd_nonrec_negacyclic.h"
cplx_ptr vector_x,vector_y,vector_res;

/******************************************************************
*
* SUPPORT CODE
*
******************************************************************/
void print_cplx(const cplx_ptr *x,int N){
  for (int i = 0; i < N; ++i)
    printf("cplxpoly[%d] = %f + i * %f\n",i,x->real[i],x->imag[i]);
  printf("\n");
}

void vec_init(){
  //PRECOMP TABLES FOR VECTOR FFT
  init_table_vctr();
  //PRECOMP FFTW
  FFTsetup();
  init_vctr();
  two_layer_init_negacyc();
  simd_init_negacyc();
  three_layer_init();
  simd_nonrec_init_negacyc();

  posix_memalign((void**)&vector_x.real,32, CPLXDIM * sizeof(double));
  posix_memalign((void**)&vector_x.imag,32, CPLXDIM * sizeof(double));
  posix_memalign((void**)&vector_y.real,32, CPLXDIM * sizeof(double));
  posix_memalign((void**)&vector_y.imag,32, CPLXDIM * sizeof(double));
  posix_memalign((void**)&vector_res.real,32, CPLXDIM * sizeof(double));
  posix_memalign((void**)&vector_res.imag,32, CPLXDIM * sizeof(double));

}

/******************************************************************
*
* FFTW MULTIPLICATION
*
******************************************************************/
void fftw_mul(ring_t *r, const ring_t *x, const ring_t *y){
  double complex cplx_x[CPLXDIM+1];
  double complex cplx_y[CPLXDIM+1];
  double complex cplx_res[CPLXDIM+1];

  FFTWforward(cplx_x,x);
  FFTWforward(cplx_y,y);

  for (int i = 0; i < CPLXDIM+1; ++i)
  {
    cplx_res[i] = cplx_x[i] * cplx_y[i];
  }
  FFTWbackward(r,cplx_res);
}

/******************************************************************
*
* FFTW DAN's TRICK MULTIPLICATION
*
******************************************************************/
void fftw_nega_mul(ring_t *r, const ring_t *x, const ring_t *y){
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  FFTW_nega_forward(cplx_x,x);
  FFTW_nega_forward(cplx_y,y);

  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }
  FFTW_nega_backward(r,cplx_res);
}

/******************************************************************
*
* NEGACYCLIC FFT RECUSIVE VECTORISED
*
******************************************************************/
void simd_negacyclic(ring_t *r, const ring_t *x, const ring_t *y)
{
  simd_forward(&vector_x,x);
  simd_forward(&vector_y,y);

  __m256d real_x,imag_x,real_y,imag_y,imag_temp,real_temp,dim;
  dim = _mm256_set1_pd(CPLXDIM);
  // double a,b,c,d;
  for (int i = 0; i < CPLXDIM; i+=4)
  {
    real_x = _mm256_load_pd(vector_x.real+i);
    imag_x = _mm256_load_pd(vector_x.imag+i);
    real_y = _mm256_load_pd(vector_y.real+i);
    imag_y = _mm256_load_pd(vector_y.imag+i);

    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    //real_temp = bd
    real_temp = _mm256_mul_pd(imag_x,imag_y);
    //imag_temp = ad
    imag_temp = _mm256_mul_pd(real_x,imag_y);
 
    real_x = _mm256_fmsub_pd(real_x,real_y,real_temp);
    imag_x = _mm256_fmadd_pd(imag_x,real_y,imag_temp);

    
    real_x = _mm256_div_pd(real_x,dim);
    imag_x = _mm256_div_pd(imag_x,dim);

    _mm256_store_pd(vector_res.real+i,real_x);
    _mm256_store_pd(vector_res.imag+i,imag_x);
  }
  simd_backward(&vector_res,r);
  // print_cplx(&vector_res,CPLXDIM);
}


/******************************************************************
*
* NEGACYCLIC FFT NONRECUSIVE VECTORISED
*
******************************************************************/
void simd_nonrec_negacyclic(ring_t *r, const ring_t *x, const ring_t *y)
{
  simd_nonrec_forward(&vector_x,x);
  simd_nonrec_forward(&vector_y,y);

  __m256d real_x,imag_x,real_y,imag_y,imag_temp,real_temp,dim;
  dim = _mm256_set1_pd(CPLXDIM);
  // double a,b,c,d;
  for (int i = 0; i < CPLXDIM; i+=4)
  {
    real_x = _mm256_load_pd(vector_x.real+i);
    imag_x = _mm256_load_pd(vector_x.imag+i);
    real_y = _mm256_load_pd(vector_y.real+i);
    imag_y = _mm256_load_pd(vector_y.imag+i);

    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    //real_temp = bd
    real_temp = _mm256_mul_pd(imag_x,imag_y);
    //imag_temp = ad
    imag_temp = _mm256_mul_pd(real_x,imag_y);
 
    real_x = _mm256_fmsub_pd(real_x,real_y,real_temp);
    imag_x = _mm256_fmadd_pd(imag_x,real_y,imag_temp);

    
    real_x = _mm256_div_pd(real_x,dim);
    imag_x = _mm256_div_pd(imag_x,dim);

    _mm256_store_pd(vector_res.real+i,real_x);
    _mm256_store_pd(vector_res.imag+i,imag_x);
  }
  simd_nonrec_backward(&vector_res,r);
  // print_cplx(&vector_res,CPLXDIM);
}

/******************************************************************
*
* NEGACYCLIC FFT MERGED LAYERS
*
******************************************************************/
void two_layer_negacyclic(ring_t *r, const ring_t *x, const ring_t *y)
{
  two_layer_forward(&vector_x,x);
  two_layer_forward(&vector_y,y);

  __m256d real_x,imag_x,real_y,imag_y,imag_temp,real_temp,dim;
  dim = _mm256_set1_pd(CPLXDIM);
  // double a,b,c,d;
  for (int i = 0; i < CPLXDIM; i+=4)
  {
    real_x = _mm256_load_pd(vector_x.real+i);
    imag_x = _mm256_load_pd(vector_x.imag+i);
    real_y = _mm256_load_pd(vector_y.real+i);
    imag_y = _mm256_load_pd(vector_y.imag+i);

    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    //real_temp = bd
    real_temp = _mm256_mul_pd(imag_x,imag_y);
    //imag_temp = ad
    imag_temp = _mm256_mul_pd(real_x,imag_y);
 
    real_x = _mm256_fmsub_pd(real_x,real_y,real_temp);
    imag_x = _mm256_fmadd_pd(imag_x,real_y,imag_temp);

    
    real_x = _mm256_div_pd(real_x,dim);
    imag_x = _mm256_div_pd(imag_x,dim);

    _mm256_store_pd(vector_res.real+i,real_x);
    _mm256_store_pd(vector_res.imag+i,imag_x);
  }
  two_layer_backward(&vector_res,r);
  // print_cplx(&vector_res,CPLXDIM);
}

/******************************************************************
*
* NEGACYCLIC FFT 3 MERGED LAYERS
*
******************************************************************/
void three_layer_negacyclic(ring_t *r, const ring_t *x, const ring_t *y)
{
  three_layer_forward(&vector_x,x);
  three_layer_forward(&vector_y,y);

  __m256d real_x,imag_x,real_y,imag_y,imag_temp,real_temp,dim;
  dim = _mm256_set1_pd(CPLXDIM);
  // double a,b,c,d;
  for (int i = 0; i < CPLXDIM; i+=4)
  {
    real_x = _mm256_load_pd(vector_x.real+i);
    imag_x = _mm256_load_pd(vector_x.imag+i);
    real_y = _mm256_load_pd(vector_y.real+i);
    imag_y = _mm256_load_pd(vector_y.imag+i);

    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    //real_temp = bd
    real_temp = _mm256_mul_pd(imag_x,imag_y);
    //imag_temp = ad
    imag_temp = _mm256_mul_pd(real_x,imag_y);
 
    real_x = _mm256_fmsub_pd(real_x,real_y,real_temp);
    imag_x = _mm256_fmadd_pd(imag_x,real_y,imag_temp);

    
    real_x = _mm256_div_pd(real_x,dim);
    imag_x = _mm256_div_pd(imag_x,dim);

    _mm256_store_pd(vector_res.real+i,real_x);
    _mm256_store_pd(vector_res.imag+i,imag_x);
  }
  three_layer_backward(&vector_res,r);
  // print_cplx(&vector_res,CPLXDIM);
}

/******************************************************************
*
* SPLIT RADIX PRECOMPUTED AND VECTORIZED NON RECURSIVE FFT MULTIPLICATION
*
******************************************************************/
void sr_vector_nonrec_mul(ring_t *r, const ring_t *x, const ring_t *y){
  fft_vector_nonrec_forward(&vector_x,x);
  fft_vector_nonrec_forward(&vector_y,y);
  __m256d real_x,imag_x,real_y,imag_y,imag_temp,real_temp;
  // double a,b,c,d;
  for (int i = 0; i < CPLXDIM; i+=4)
  {
    real_x = _mm256_load_pd(vector_x.real+i);
    imag_x = _mm256_load_pd(vector_x.imag+i);
    real_y = _mm256_load_pd(vector_y.real+i);
    imag_y = _mm256_load_pd(vector_y.imag+i);

    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    //real_temp = bd
    real_temp = _mm256_mul_pd(imag_x,imag_y);
    //imag_temp = ad
    imag_temp = _mm256_mul_pd(real_x,imag_y);

    real_x = _mm256_fmsub_pd(real_x,real_y,real_temp);
    imag_x = _mm256_fmadd_pd(imag_x,real_y,imag_temp);

    real_y = _mm256_set1_pd(CPLXDIM);
    real_x = _mm256_div_pd(real_x,real_y);
    imag_x = _mm256_div_pd(imag_x,real_y);

    _mm256_store_pd(vector_res.real+i,real_x);
    _mm256_store_pd(vector_res.imag+i,imag_x);

  }
  fft_vector_nonrec_backward(&vector_res,r);
}
/******************************************************************
*
* SPLIT RADIX PRECOMPUTED AND VECTORIZED FFT MULTIPLICATION
*
******************************************************************/
void sr_vector_mul(ring_t *r, const ring_t *x, const ring_t *y){
  // printf("\n\n**************split-radix FAST**************\n");

  fft_vector_forward(&vector_x,x);
  fft_vector_forward(&vector_y,y);
  
  __m256d real_x,imag_x,real_y,imag_y,imag_temp,real_temp;
  // double a,b,c,d;
  for (int i = 0; i < CPLXDIM; i+=4)
  {
    real_x = _mm256_load_pd(vector_x.real+i);
    imag_x = _mm256_load_pd(vector_x.imag+i);
    real_y = _mm256_load_pd(vector_y.real+i);
    imag_y = _mm256_load_pd(vector_y.imag+i);

    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    //real_temp = bd
    real_temp = _mm256_mul_pd(imag_x,imag_y);
    //imag_temp = ad
    imag_temp = _mm256_mul_pd(real_x,imag_y);
     
    real_x = _mm256_fmsub_pd(real_x,real_y,real_temp);
    imag_x = _mm256_fmadd_pd(imag_x,real_y,imag_temp);

    real_y = _mm256_set1_pd(CPLXDIM);
    real_x = _mm256_div_pd(real_x,real_y);
    imag_x = _mm256_div_pd(imag_x,real_y);

    _mm256_store_pd(vector_res.real+i,real_x);
    _mm256_store_pd(vector_res.imag+i,imag_x);
  }
  fft_vector_backward(&vector_res,r);
}
