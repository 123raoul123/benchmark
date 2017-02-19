#ifndef LUT_NEGACYCLIC_H
#define LUT_NEGACYCLIC_H

#include <complex.h>

void init_negacyc();
void inverse_phi_lut(double complex *x,int n,int lo,int level);
void recursive_phi_lut(double complex *x,int n,int lo,int level);

#endif