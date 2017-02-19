#ifndef TEST_H
#define TEST_H

#define CPLXDIM 512
#define REALDIM (2*CPLXDIM)
#define ZEROPAD (2*REALDIM)
#define ROOTDIM (2*REALDIM)

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define W(N,k) (cexp(2.0 * M_PI * I * (double)k / (double) N))

#include <stdint.h>
#include <complex.h>

typedef struct {
	double real[CPLXDIM];
	double imag[CPLXDIM];
} cplx;

typedef struct {
  uint32_t v[REALDIM];
} ring_t;

#endif