#ifndef VEC_MUL_H
#define VEC_MUL_H

#define W(N,k) (cexp(2.0 * M_PI * I * (double)k / (double) N))
#define calc_cos(N,k) (cos(2.0 * M_PI * (double)k / (double) N))
#define calc_sin(N,k) (sin(2.0 * M_PI * (double)k / (double) N))

typedef struct
{	
	double *real;
	double *imag;
}cplx_ptr;

void vec_init();
void fftw_nega_mul(ring_t *r, const ring_t *x, const ring_t *y);
void fftw_mul(ring_t *r, const ring_t *x, const ring_t *y);
void naive_real_mul(ring_t *r, const ring_t *x, const ring_t *y);
void sr_vector_mul(ring_t *r, const ring_t *x, const ring_t *y);
void sr_vector_nonrec_mul(ring_t *r, const ring_t *x, const ring_t *y);
void two_layer_negacyclic(ring_t *r, const ring_t *x, const ring_t *y);
void simd_negacyclic(ring_t *r, const ring_t *x, const ring_t *y);
void three_layer_negacyclic(ring_t *r, const ring_t *x, const ring_t *y);
void simd_nonrec_negacyclic(ring_t *r, const ring_t *x, const ring_t *y);
void simd_nonmerge(ring_t *r, const ring_t *x, const ring_t *y);
void nonrec_nonmerge(ring_t *r, const ring_t *x, const ring_t *y);

#endif