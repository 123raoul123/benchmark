#ifndef LUT_SPLIT_RADIX_H
#define LUT_SPLIT_RADIX_H

void fft_precompsr_forward(cplx *x);
void fft_precompsr_backward(cplx *x,ring_t *r);
void init_table();
void destruct_table();

#endif