#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <immintrin.h>
#include "../../test.h"
#include "../vec_mul.h"
#include "simd_experimental_negacyclic.h"


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


void simd_experimental_init_negacyc()
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
static void inverse_simd(cplx_ptr *x)
{	
	__m256d real_x,imag_x,real_y,imag_y,imag_twid,real_twid,temp_real,temp_imag,sub_real,sub_imag;
	int m = 4;
	int step = 8;
	int level;

	for(int lo =0;lo < CPLXDIM;lo = lo+step)
	{
		level = 8;
		real_x = _mm256_load_pd(x->real+lo);
		imag_x = _mm256_load_pd(x->imag+lo);
		real_y = _mm256_load_pd(x->real+lo+m);
		imag_y = _mm256_load_pd(x->imag+lo+m);
		//LAYER 2
		//we need to collect all the left and the right parts
		//real part
		sub_real = _mm256_unpacklo_pd(real_x,real_y);
		sub_imag = _mm256_unpackhi_pd(real_x,real_y);

		real_y = _mm256_add_pd(sub_real,sub_imag);
		real_x = _mm256_sub_pd(sub_real,sub_imag);

		//imag part
		sub_real = _mm256_unpacklo_pd(imag_x,imag_y);
		sub_imag = _mm256_unpackhi_pd(imag_x,imag_y);

		imag_y = _mm256_add_pd(sub_real,sub_imag);
		imag_x = _mm256_sub_pd(sub_real,sub_imag);

		//Now do complex inverse computation with roots of unity
		real_twid = _mm256_setr_pd(vec_wortel[0][level][lo/2],vec_wortel[0][level][(lo+m)/2],vec_wortel[0][level][(lo+2)/2],vec_wortel[0][level][(lo+m+2)/2]);
    	imag_twid = _mm256_setr_pd(vec_wortel[2][level][lo/2],vec_wortel[2][level][(lo+m)/2],vec_wortel[2][level][(lo+2)/2],vec_wortel[2][level][(lo+m+2)/2]);
	    temp_real = _mm256_mul_pd(imag_x,imag_twid);
    	temp_imag = _mm256_mul_pd(imag_x,real_twid);

		//TEMP_imag = ad + bc
		imag_x = _mm256_fmadd_pd(real_x,imag_twid,temp_imag);
		//TEMP_real = ac - bd
		real_x = _mm256_fmsub_pd(real_x,real_twid,temp_real);

		//LAYER 4
		--level;
		//real part
		sub_real = _mm256_permute2f128_pd(real_y,real_x,0x20);
		sub_imag = _mm256_permute2f128_pd(real_y,real_x,0x31);

		real_y = _mm256_add_pd(sub_real,sub_imag);
		real_x = _mm256_sub_pd(sub_real,sub_imag);
		//imag part
		sub_real = _mm256_permute2f128_pd(imag_y,imag_x,0x20);
		sub_imag = _mm256_permute2f128_pd(imag_y,imag_x,0x31);

		imag_y = _mm256_add_pd(sub_real,sub_imag);
		imag_x = _mm256_sub_pd(sub_real,sub_imag);
		//mult
		real_twid = _mm256_setr_pd(vec_wortel[0][level][lo/4],vec_wortel[0][level][(lo+m)/4],vec_wortel[0][level][lo/4],vec_wortel[0][level][(lo+m)/4]);
	    imag_twid = _mm256_setr_pd(vec_wortel[2][level][lo/4],vec_wortel[2][level][(lo+m)/4],vec_wortel[2][level][lo/4],vec_wortel[2][level][(lo+m)/4]);
	    temp_real = _mm256_mul_pd(imag_x,imag_twid);
    	temp_imag = _mm256_mul_pd(imag_x,real_twid);
		//TEMP_imag = ad + bc
		imag_x = _mm256_fmadd_pd(real_x,imag_twid,temp_imag);
		//TEMP_real = ac - bd
		real_x = _mm256_fmsub_pd(real_x,real_twid,temp_real);
		
		//LAYER 8
		--level;
		//real part
		temp_real = _mm256_unpacklo_pd(real_y,real_x);
		temp_imag = _mm256_unpackhi_pd(real_y,real_x);

		real_x = _mm256_add_pd(temp_real,temp_imag);
		real_y = _mm256_sub_pd(temp_real,temp_imag);
		//imag part
		sub_real  = _mm256_unpacklo_pd(imag_y,imag_x);
		sub_imag  = _mm256_unpackhi_pd(imag_y,imag_x);

		imag_x = _mm256_add_pd(sub_real,sub_imag);
		imag_y = _mm256_sub_pd(sub_real,sub_imag);
		//mult
	    real_twid = _mm256_set1_pd(vec_wortel[0][level][lo/step]);
    	imag_twid = _mm256_set1_pd(vec_wortel[2][level][lo/step]);
	    temp_real = _mm256_mul_pd(imag_y,imag_twid);
    	temp_imag = _mm256_mul_pd(imag_y,real_twid);

		//TEMP_imag = ad + bc
		imag_y = _mm256_fmadd_pd(real_y,imag_twid,temp_imag);
		//TEMP_real = ac - bd
		real_y = _mm256_fmsub_pd(real_y,real_twid,temp_real);

		real_x =_mm256_permute4x64_pd(real_x,0xd8);
		real_y =_mm256_permute4x64_pd(real_y,0xd8);
		imag_x =_mm256_permute4x64_pd(imag_x,0xd8);
		imag_y =_mm256_permute4x64_pd(imag_y,0xd8);

		_mm256_store_pd(x->real+lo,real_x);
		_mm256_store_pd(x->imag+lo,imag_x);
		_mm256_store_pd(x->real+lo+m,real_y);
		_mm256_store_pd(x->imag+lo+m,imag_y);
	}

	m = 8;
	step = 16;
  	for(level = 5;level >= 0; --level)
  	{
	  	for(int lo =0; lo<CPLXDIM;lo = lo+step)
		{

			real_twid = _mm256_set1_pd(vec_wortel[0][level][lo/step]);
			imag_twid = _mm256_set1_pd(vec_wortel[2][level][lo/step]);
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
		m = step;
		step = step*2;
	}
}

static void iterative_simd(cplx_ptr *x)
{
	__m256d real_x,imag_x,real_y,imag_y,imag_twid,real_twid,temp_real,temp_imag;
  	int m = CPLXDIM/2;
  	int step = CPLXDIM;
  	int level = 0;
  	for(;level < 7; ++level)
  	{
	  	for(int lo =0; lo<CPLXDIM;lo = lo+step)
		{
			real_twid = _mm256_set1_pd(vec_wortel[0][level][lo/step]);
			imag_twid = _mm256_set1_pd(vec_wortel[1][level][lo/step]);

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
		}
		step = m;
		m = m/2;
  	}
  	__m128d omg_real,omg_imag,x_real,x_imag,y_real,y_imag,real_temp,imag_temp;
  	for(int lo = 0; lo<CPLXDIM;lo= lo+step)
  	{	
  		x_real = _mm_load_pd(x->real+lo);
  		x_imag = _mm_load_pd(x->imag+lo);
  		y_real = _mm_load_pd(x->real+lo+m);
  		y_imag = _mm_load_pd(x->imag+lo+m);
  		omg_real = _mm_set1_pd(vec_wortel[0][level][lo/4]);
  		omg_imag = _mm_set1_pd(vec_wortel[1][level][lo/4]);

  	
		real_temp = _mm_mul_pd(y_imag,omg_imag);
		imag_temp = _mm_mul_pd(y_imag,omg_real);
		//TEMP_real = ac - bd
		real_temp = _mm_fmsub_pd(y_real,omg_real,real_temp);
		//TEMP_imag = ad + bc
		imag_temp = _mm_fmadd_pd(y_real,omg_imag,imag_temp);

		y_real = _mm_sub_pd(x_real,real_temp);
		y_imag = _mm_sub_pd(x_imag,imag_temp);

		x_real = _mm_add_pd(x_real,real_temp);
		x_imag = _mm_add_pd(x_imag,imag_temp);

		_mm_store_pd(x->real+lo,x_real);
		_mm_store_pd(x->imag+lo,x_imag);
		_mm_store_pd(x->real+lo+m,y_real);
		_mm_store_pd(x->imag+lo+m,y_imag);
	}
	++level;
	double i_twid,r_twid,i_temp,r_temp;
	for(int lo=0;lo<CPLXDIM;lo=lo+2)
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

void simd_experimental_forward(cplx_ptr *x,const ring_t *ring)
{
  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    x->real[i] = ring->v[i];
    x->imag[i] = ring->v[j];
    ++j;
  }
  iterative_simd(x);
}

void simd_experimental_backward(cplx_ptr *x, ring_t *ring)
{
  inverse_simd(x);

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    ring->v[i] = x->real[i];
    ring->v[j] = x->imag[i];
    ++j;
  }
}