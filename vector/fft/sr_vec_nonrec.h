#ifndef SR_VEC_NONREC_H
#define SR_VEC_NONREC_H

void init_vctr();
void destruct_vctr();
void fft_vector_nonrec_forward(cplx_ptr *x,const ring_t *ring);
void fft_vector_nonrec_backward(cplx_ptr *cplx_x,ring_t *res);

#endif
