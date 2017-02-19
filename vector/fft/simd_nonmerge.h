#ifndef SIMD_NONMERGE_H
#define SIMD_NONMERGE_H

void nonmerge_init();
void nonmerge_forward(cplx_ptr *x,const ring_t *ring);
void nonmerge_backward(cplx_ptr *x, ring_t *ring);

#endif