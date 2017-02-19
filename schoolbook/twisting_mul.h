#ifndef TWISTING_MUL_H
#define TWISTING_MUL_H

void twist_tangent(ring_t *r, const ring_t *x, const ring_t *y);
void twist_split_radix(ring_t *r, const ring_t *x, const ring_t *y);
void twist_twisted(ring_t *r, const ring_t *x, const ring_t *y);
void twist_negacyclic(ring_t *r, const ring_t *x, const ring_t *y);
void twist_fiduccia(ring_t *r, const ring_t *x, const ring_t *y);
#endif