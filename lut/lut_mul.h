#ifndef MUL_H
#define MUL_H

void lut_init();
void lut_tangent(ring_t *r, const ring_t *x, const ring_t *y);
void lut_split_radix(ring_t *r, const ring_t *x, const ring_t *y);
void lut_negacyclic(ring_t *r, const ring_t *x, const ring_t *y);
#endif