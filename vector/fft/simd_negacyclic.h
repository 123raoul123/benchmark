#ifndef SIMD_NEGACYCLIC_H
#define SIMD_NEGACYCLIC_H

void simd_init_negacyc();
void simd_forward(cplx_ptr *x,const ring_t *ring);
void simd_backward(cplx_ptr *x, ring_t *ring);

#endif