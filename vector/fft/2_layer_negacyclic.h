#ifndef TWO_LAYER_NEGACYCLIC_H
#define TWO_LAYER_NEGACYCLIC_H

void two_layer_init_negacyc();
void two_layer_forward(cplx_ptr *x,const ring_t *ring);
void two_layer_backward(cplx_ptr *x, ring_t *ring);

#endif