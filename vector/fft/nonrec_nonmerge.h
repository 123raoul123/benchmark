#ifndef NONREC_NONMERGE_H
#define NONREC_NONMERGE_H

void nonrec_nonmerge_init();
void nonrec_nonmerge_forward(cplx_ptr *x,const ring_t *ring);
void nonrec_nonmerge_backward(cplx_ptr *x, ring_t *ring);

#endif