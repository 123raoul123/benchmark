#ifndef SIMD_EXPERIMENTAL_NEGACYCLIC_H
#define SIMD_EXPERIMENTAL_NEGACYCLIC_H

void simd_experimental_init_negacyc();
void simd_experimental_forward(cplx_ptr *x,const ring_t *ring);
void simd_experimental_backward(cplx_ptr *x, ring_t *ring);

#endif