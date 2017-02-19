#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <immintrin.h>
#include "../../test.h"
#include "../vec_mul.h"
#include "simd_negacyclic.h"


static double ***vec_wortel;

static void init_simd_wortel(int n,int lo,int level,double complex root)
{	
	if(n > 1){
		int m = n/2;
		vec_wortel[0][level][lo/n] = creal(root);
		vec_wortel[1][level][lo/n] = cimag(root);
		vec_wortel[2][level][lo/n] = -cimag(root);
		++level;
		init_simd_wortel(m,lo,level,csqrt(root));
		init_simd_wortel(m,lo+m,level,csqrt(-root));
	}
}


void nonmerge_init()
{	
	vec_wortel = malloc(3*sizeof(**vec_wortel));
	int loga = log2(CPLXDIM);
	vec_wortel[0] = malloc(loga*sizeof(*vec_wortel));
	vec_wortel[1] = malloc(loga*sizeof(*vec_wortel));
	vec_wortel[2] = malloc(loga*sizeof(*vec_wortel));
	int j = 1;
	for (int i = 0; i < loga; ++i)
	{
		vec_wortel[0][i] = malloc(j*sizeof(double));
		vec_wortel[1][i] = malloc(j*sizeof(double));
		vec_wortel[2][i] = malloc(j*sizeof(double));
		j = j <<1;
	}
	init_simd_wortel(CPLXDIM,0,0,csqrt(I));
}

/******************************************************************
*
* SMART COMPLEX MULTIPLICATION
*
******************************************************************/
static void inverse_simd(cplx_ptr *x,int n,int lo,int level)
{	
	int m = n/2;
	if(n > 4)
	{
		inverse_simd(x,m,lo,level+1);
		inverse_simd(x,m,lo+m,level+1);
	    __m256d real_x,imag_x,real_y,imag_y,imag_twid,real_twid,temp_real,temp_imag;
    	real_twid = _mm256_set1_pd(vec_wortel[0][level][lo/n]);
    	imag_twid = _mm256_set1_pd(vec_wortel[2][level][lo/n]);
		for(int i=lo;i<m+lo;i+=4)
		{	
			real_x = _mm256_load_pd(x->real+i);
		    imag_x = _mm256_load_pd(x->imag+i);
		    real_y = _mm256_load_pd(x->real+i+m);
		    imag_y = _mm256_load_pd(x->imag+i+m);

		    temp_real = real_x;
		    temp_imag = imag_x;

		    real_x = _mm256_add_pd(temp_real,real_y);
	  		imag_x = _mm256_add_pd(temp_imag,imag_y);

	  		real_y = _mm256_sub_pd(temp_real,real_y);
	  		imag_y = _mm256_sub_pd(temp_imag,imag_y);

	  		//TEMP_real = bd
		    temp_real = _mm256_mul_pd(imag_y,imag_twid);
		    //TEMP_imag = bc
		    temp_imag = _mm256_mul_pd(imag_y,real_twid);

		    //imag_y = ad + bc
			imag_y = _mm256_fmadd_pd(real_y,imag_twid,temp_imag);
		    //real_y = ac - bd
		    real_y = _mm256_fmsub_pd(real_y,real_twid,temp_real);
			

			_mm256_store_pd(x->real+i,real_x);
		    _mm256_store_pd(x->imag+i,imag_x);
		    _mm256_store_pd(x->real+i+m,real_y);
		    _mm256_store_pd(x->imag+i+m,imag_y);
		}
	}
	double r_twid,i_twid,r_temp,i_temp,a,b;
	if(n == 4)
	{	
		inverse_simd(x,m,lo,level+1);
		inverse_simd(x,m,lo+m,level+1);

		r_twid = vec_wortel[0][level][lo/n];
		i_twid = vec_wortel[2][level][lo/n];

		for(int i=lo;i<lo+m;++i)
		{
			a = x->real[i]-x->real[i+m];
			b = x->imag[i]-x->imag[i+m];

			r_temp = (a*r_twid) - (b*i_twid);
			i_temp = (a*i_twid) + (b*r_twid);

			x->real[i] = x->real[i]+x->real[i+m];
			x->imag[i] = x->imag[i]+x->imag[i+m];

			x->real[i+m] = r_temp;
			x->imag[i+m] = i_temp;
		}

	}
	if(n == 2)
	{
		r_twid = vec_wortel[0][level][lo/n];
		i_twid = vec_wortel[2][level][lo/n];

		a = x->real[lo]-x->real[lo+m];
		b = x->imag[lo]-x->imag[lo+m];

		r_temp = (a*r_twid) - (b*i_twid);
		i_temp = (a*i_twid) + (b*r_twid);

		x->real[lo] = x->real[lo]+x->real[lo+m];
		x->imag[lo] = x->imag[lo]+x->imag[lo+m];

		x->real[lo+m] = r_temp;
		x->imag[lo+m] = i_temp;
	}
}

static void recursive_simd(cplx_ptr *x,int n,int lo,int level)
{	
  int m = n/2;	
  if(n > 4)
  {
    __m256d real_x,imag_x,real_y,imag_y,imag_twid,real_twid,temp_real,temp_imag;
    real_twid = _mm256_set1_pd(vec_wortel[0][level][lo/n]);
    imag_twid = _mm256_set1_pd(vec_wortel[1][level][lo/n]);
    for(int i=lo; i < lo+m;i+=4)
    {
	  //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
	  real_x = _mm256_load_pd(x->real+i);
      imag_x = _mm256_load_pd(x->imag+i);
      real_y = _mm256_load_pd(x->real+i+m);
      imag_y = _mm256_load_pd(x->imag+i+m);
      //TEMP_real = bd
      temp_real = _mm256_mul_pd(imag_y,imag_twid);
      //TEMP_imag = bc
      temp_imag = _mm256_mul_pd(imag_y,real_twid);

      //TEMP_real = ac - bd
      temp_real = _mm256_fmsub_pd(real_y,real_twid,temp_real);
	  //TEMP_imag = ad + bc
	  temp_imag = _mm256_fmadd_pd(real_y,imag_twid,temp_imag);

	  real_y = _mm256_sub_pd(real_x,temp_real);
	  imag_y = _mm256_sub_pd(imag_x,temp_imag);

	  real_x = _mm256_add_pd(real_x,temp_real);
	  imag_x = _mm256_add_pd(imag_x,temp_imag);

	  _mm256_store_pd(x->real+i,real_x);
      _mm256_store_pd(x->imag+i,imag_x);
      _mm256_store_pd(x->real+i+m,real_y);
      _mm256_store_pd(x->imag+i+m,imag_y);
	}
	++level;
    recursive_simd(x,m,lo,level);
    recursive_simd(x,m,lo + m,level);
  }
  double r_twid,i_twid,r_temp,i_temp;
  if(n == 4)
  {
	r_twid = vec_wortel[0][level][lo/n];
	i_twid = vec_wortel[1][level][lo/n];

  	for(int i=lo;i<lo+m;++i)
  	{
		r_temp = x->imag[i+m] * i_twid;
		i_temp = x->imag[i+m] * r_twid;

		r_temp = (x->real[i+m] * r_twid)-r_temp;
		i_temp = (x->real[i+m] * i_twid)+i_temp;

		x->real[i+m] = x->real[i] - r_temp;
		x->imag[i+m] = x->imag[i] - i_temp;

		x->real[i] = x->real[i] + r_temp;
		x->imag[i] = x->imag[i] + i_temp;

  	}
  	++level;
    recursive_simd(x,m,lo,level);
    recursive_simd(x,m,lo + m,level);
  }
  if(n == 2)
  {
	r_twid = vec_wortel[0][level][lo/2];
	i_twid = vec_wortel[1][level][lo/2];

	r_temp = x->imag[lo+1] * i_twid;
	i_temp = x->imag[lo+1] * r_twid;

	r_temp = (x->real[lo+1] * r_twid)-r_temp;
	i_temp = (x->real[lo+1] * i_twid)+i_temp;

	x->real[lo+1] = x->real[lo] - r_temp;
	x->imag[lo+1] = x->imag[lo] - i_temp;

	x->real[lo] = x->real[lo] + r_temp;
	x->imag[lo] = x->imag[lo] + i_temp;
  }
}

void nonmerge_forward(cplx_ptr *x,const ring_t *ring)
{
  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    x->real[i] = ring->v[i];
    x->imag[i] = ring->v[j];
    ++j;
  }
  recursive_simd(x, CPLXDIM,0,0);
}

void nonmerge_backward(cplx_ptr *x, ring_t *ring)
{
  inverse_simd(x, CPLXDIM,0,0);

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    ring->v[i] = x->real[i];
    ring->v[j] = x->imag[i];
    ++j;
  }
}