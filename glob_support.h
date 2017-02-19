#ifndef GLOB_SUPPORT_H
#define GLOB_SUPPORT_H

void to_real(const double complex *cplx_x, ring_t *x);
void to_complex(const ring_t *x, double complex *cplx_x);
void print_complex(const double complex *a, int N);
void print_double(const cplx *x,int N);

#endif