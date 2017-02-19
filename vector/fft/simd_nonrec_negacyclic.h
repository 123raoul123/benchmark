#ifndef SIMD_NONREC_NEGACYCLIC_H
#define SIMD_NONREC_NEGACYCLIC_H

void simd_nonrec_init_negacyc();
void simd_nonrec_forward(cplx_ptr *x,const ring_t *ring);
void simd_nonrec_backward(cplx_ptr *x, ring_t *ring);

#endif