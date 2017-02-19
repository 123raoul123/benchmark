#ifndef ZEROPAD_MUL_H
#define ZEROPAD_MUL_H

void zeropad_tangent(ring_t *r, const ring_t *x, const ring_t *y);
void zeropad_split_radix(ring_t *r, const ring_t *x, const ring_t *y);
void zeropad_fiduccia(ring_t *r, const ring_t *x, const ring_t *y);
void zeropad_twisted(ring_t *r, const ring_t *x, const ring_t *y);
void schoolbook(ring_t *r, const ring_t *x, const ring_t *y);
#endif