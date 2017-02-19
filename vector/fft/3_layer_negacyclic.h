#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <immintrin.h>
#include "../vec_mul.h"


void three_layer_forward(cplx_ptr *x,const ring_t *ring);
void three_layer_backward(cplx_ptr *x, ring_t *ring);
void three_layer_init();