#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lut_tangent.h"
#include "../../test.h"


double complex **twidle;
double **base;


double lut_s(int n, int k) 
{
  if (n <= 4) 
    return 1.0;

  int k4 = k % (n/4);

  if (k4 <= n/8) 
    return (lut_s(n/4,k4) * cos(2.0 * M_PI * (double)k4 / (double)n));
  
  return (lut_s(n/4,k4) * sin(2.0f * M_PI * (double)k4 / (double)n));
}

void init_tangent()
{
  base = malloc(11*sizeof(*base));
  int slide = 1;
  for (int i = 1; i < 11; ++i)
  {
    // printf("i = %d \n", );
    base[i] = malloc(CPLXDIM*sizeof(double));
    for (int j = 0; j < CPLXDIM; ++j)
    {
      base[i][j] = lut_s(slide,j);
    }
    slide = slide << 1;
  }
  twidle = malloc(2*sizeof(*twidle));
  twidle[0] = malloc(ROOTDIM*sizeof(double complex));
  twidle[1] = malloc(ROOTDIM*sizeof(double complex));
  for (int i = 0; i < ROOTDIM; ++i)
  {
    twidle[0][i] = W(ROOTDIM,i);
    twidle[1][i] = conj(twidle[0][i]);
  }
}

void lut_twist(double complex *cplx_x,int n,int m,int lo)
{
  // printf("n = %d, m = %d, lo = %d\n",n,m,lo );
  int j = 1,scale = ROOTDIM/n;
  for (int i = lo+1; i < lo+m; ++i)
  {
    // printf("i = %d, j = %d\n",i,j);
    cplx_x[i] = cplx_x[i] * twidle[0][j*scale];
    ++j;
  }
}

void lut_untwist(double complex *cplx_x,int n,int m,int lo)
{
  // printf("n = %d, m = %d, lo = %d\n",n,m,lo );
  int j = 1,scale = ROOTDIM/n;
  for (int i = lo+1; i < lo+m; ++i)
  {
    // printf("i = %d, j = %d\n",i,j);
    cplx_x[i] = cplx_x[i] * twidle[1][j*scale];
    ++j;
  }
}

void lut_base_change(double complex *x,int n,int n_new,int m,int lo)
{ int j = 1,scale = log2(n)+1;
  if(n_new == -1)
  {
    for (int i = lo+1; i < lo+m; ++i)
    {
      x[i] = x[i] * base[scale][j];
      ++j;
    }
  }
  else
  {
    int new_scale = log2(n_new)+1;
    double temp;
    for (int i = lo+1; i < lo+m; ++i)
    { 
      temp = base[new_scale][j]/base[scale][j];
      // printf("temp = %f, x[%d] = %f + I %f\n",temp,i,creal(x[i]),cimag(x[i]));
      x[i] = x[i] * temp;
      ++j;
    }
  }
}

void lut_cost_4_twist(double complex *cplx_x,int n,int m,int lo)
{  
  int j = 1,scale = ROOTDIM/n;
  int new_slide = log2(m)+1,slide = log2(n)+1;
  double complex temp;
  // printf("lo = %d m = %d n = %d\n",lo,m,n);
  for (int i = lo+1; i < lo+m; ++i)
  {
    temp = base[new_slide][j]/base[slide][j];
    // printf("temp = %f + i %f\n",creal(temp),cimag(temp));
    // printf("W(%d,%d) = %f + i %f\n",n,j,creal(W(n,j)),cimag(W(n,j)));
    // printf("temp = %f + i %f, x[%d] = %f + I %f\n",creal(temp),cimag(temp),i,creal(cplx_x[i]),cimag(cplx_x[i]));
    cplx_x[i] = cplx_x[i] * twidle[0][scale*j] *temp;
    ++j;
  }
}

void lut_cost_4_untwist(double complex *cplx_x,int n,int m,int lo)
{
  int j = 1, scale = ROOTDIM/n;
  int new_slide = log2(m)+1,slide = log2(n)+1;
  double complex temp;
  for (int i = lo+1; i < lo+m; ++i)
  {
    temp = base[new_slide][j]/base[slide][j];
    cplx_x[i] = cplx_x[i] * twidle[1][scale*j]*temp;
    ++j;
  }
}

void lut_cost_4_untwist_inverse(double complex *cplx_x,int n,int m,int lo)
{
  int j = 1, scale = ROOTDIM/n;
  int new_slide = log2(m)+1,slide = log2(n)+1;
  double complex temp;
  for (int i = lo+1; i < lo+m; ++i)
  {
    temp = base[slide][j]/base[new_slide][j];
    // temp = lut_s(n,j)/lut_s(m,j);
    cplx_x[i] = cplx_x[i] *  twidle[1][scale*j]*temp;
    ++j;
  }
}

void lut_cost_4_twist_inverse(double complex *cplx_x,int n,int m,int lo)
{
  int j = 1, scale = ROOTDIM/n;
  int new_slide = log2(m)+1,slide = log2(n)+1;
  double complex temp;
  // printf("lo = %d m = %d n = %d\n",lo,m,n);
  for (int i = lo+1; i < lo+m; ++i)
  {
    temp = base[slide][j]/base[new_slide][j];
    // temp = lut_s(n,j)/lut_s(m,j);
    // printf("temp = %f + i %f\n",creal(temp),cimag(temp));
    // printf("W(%d,%d) = %f + i %f\n",n,j,creal(W(n,j)),cimag(W(n,j)));
    // printf("temp = %f + i %f, x[%d] = %f + I %f\n",creal(temp),cimag(temp),i,creal(cplx_x[i]),cimag(cplx_x[i]));
    cplx_x[i] = cplx_x[i] * twidle[0][scale*j] *temp;
    ++j;
  }
}

void lut_tangent_8(double complex *x,int n,int lo)
{ 
  double complex temp;
  if(n == 2){
    // printf("LEVEL2\n");
    temp = x[lo];
    x[lo] = temp + x[lo+1];
    x[lo+1] = temp - x[lo+1];
    // print_complex(x,8);
  }
  else if(n == 4)
  {
    int m = n/2;
    //Go from (x^4 +1 to x^2 -1 and x^2 +1)
    for(int i=lo; i < lo+m;++i){
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }
    lut_tangent_8(x,m,lo);
  lo = lo+m;

    //Go from (x^2 +1 to x -i and x +i)
  temp = x[lo];
  x[lo] = temp + I * x[lo+1];
  x[lo+1] = temp - I * x[lo+1];
  }
  else if(n>=8){
     int m = n/2;
      //Go from (x^8n +1 to x^4n -1 and x^4n +1)
      for(int i=lo; i < lo+m;++i){
        temp = x[i];
        x[i] = temp + x[i+m];
        x[i+m] = temp - x[i+m];
      }
      //NOW take the left side
      //Go from (x^4n-1) to (x^2n-1) and (x^2n+1)
      int right_lo = lo + m;
      m = m/2;
      for(int i=lo; i < lo+m;++i){
        temp = x[i];
        x[i] = temp + x[i+m];
        x[i+m] = temp - x[i+m];
      }
      //NOW take the right side
      //Go from (x^4n+1) to (x^2n-i) and (x^2n+i)
      for (int i = right_lo; i < right_lo+m; ++i)
      {
        // printf("i = %d, m = %d\n",i,m );
        temp = x[i];
        x[i] = temp + I * x[i+m];
        x[i+m] = temp - I * x[i+m];
      }
      //BEFORE CHANGING M SET OTHER LO VALUES FOR LEFT_RIGHT AND RIGHT_RIGHT BRANCH
      //DO COST-4 TWIST
      lut_cost_4_twist(x,n,m,right_lo);
      lut_tangent_8(x,m,right_lo);
      //DO COST-4 UNTWIST
      lut_cost_4_untwist(x,n,m,right_lo+m);
      lut_tangent_8(x,m,right_lo+m);
      //BACK TO THE LEFT
      //NOW TO A BASE TWIST GO FROM X^2n-1 BASE lut_s(8n,k) to X^2n-1 BASE lut_s(2n,k)
      lut_base_change(x,n,m,m,lo);
      lut_tangent_8(x,m,lo);
      lo = lo+m;
      //FINISH THE LEFT_RIGHT BRANCH
      //ALSO BASE CHANGE FOR MID GO FROM X^2n+1 BASE lut_s(8n,k) to X^2n+1 BASE lut_s(4n,k)
      // printf("BASECHANGE!!!\n");
      // print_complex(x,8);
      lut_base_change(x,n,n/2,m,lo);
      // printf("BASECHANGE!!!\n");
      // print_complex(x,8);
      m=m/2;
      //GO FROM (X^2n+1) TO (X^n-i) AND (x^n+i)
      for (int i = lo; i < lo+m; ++i)
      {
        // printf("i = %d, m = %d\n",i,m );
        temp = x[i];
        x[i] = temp + I * x[i+m];
        x[i+m] = temp - I * x[i+m];
      }
      // lut_base_change(x,lo,m,n/2,n/8);
      // lut_twist(x,n/2,m,lo);
      lut_cost_4_twist(x,n/2,m,lo);
      lut_tangent_8(x,m,lo);
      // lut_base_change(x,lo+m,m,n/2,n/8);
      // lut_untwist(x,n/2,m,lo+m);
      lut_cost_4_untwist(x,n/2,m,lo+m);
      lut_tangent_8(x,m,lo+m);
  }
}

void lut_tangent_8_inverse(double complex *x,int n,int lo)
{
  double complex temp;
  if(n == 2){
    // printf("LEVEL2\n");
    temp = x[lo];
    x[lo] = temp + x[lo+1];
    x[lo+1] = temp - x[lo+1];
    // print_complex(x,8);
  }
  else if(n == 4)
  {
    int m = n/2;
    lo = lo +m;
    //Go from  x-i and x +i to (x^2 +1)
  temp = x[lo];
  x[lo] = temp + x[lo+1];
  x[lo+1] = -I*(temp - x[lo+1]);
  lo = lo -m;
    //GET LEFT PART
    lut_tangent_8(x,m,lo);
    //Go from x^2 -1 and x^2 +1 to x^4 +1
    for(int i=lo; i < lo+m;++i){
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }
  }
  else if(n>=8){
    int right_lo = lo + n/2;
    int m = n/4;
    lo = lo+m;
    m = m/2;
    lut_tangent_8_inverse(x,m,lo+m);
    lut_cost_4_twist_inverse(x,n/2,m,lo+m);
    lut_tangent_8_inverse(x,m,lo);
    lut_cost_4_untwist_inverse(x,n/2,m,lo);
    //GO FROM (X^n-i) AND (x^n+i) TO (X^2n+1)
    for (int i = lo; i < lo+m; ++i)
    {
      // printf("i = %d, m = %d\n",i,m );
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = -I*(temp - x[i+m]);
    }
    m = m*2;
    //NOW WE NEED TO CHANGE FROM BASE lut_s(4n,k) TO lut_s(8n,k)
    // printf("BASECHANGE!!!\n");
    // print_complex(x,8);
    lut_base_change(x,n/2,n,m,lo);
    // printf("BASECHANGE!!!\n");
    // print_complex(x,8);
    lo = lo-m;
    //BACK TO THE LEFT WE FIRST NEED TO GET LEFT BRANCH (x^2n-1) BASE lut_s(2n,k)
    lut_tangent_8_inverse(x,m,lo);
    //NOW TO A BASE CHANGE GO FROM X^2n-1 BASE lut_s(2n,k) to X^2n-1 BASE lut_s(8n,k)
    lut_base_change(x,n/4,n,m,lo);
    //NOW GO BACK TO THE RIGHT PART
    //GET RIGHT (x^2n-1)
    lut_tangent_8_inverse(x,m,right_lo+m);
    //TWIST (x^2n-1) to (x^2n+i) with base lut_s(8n,k)
    lut_cost_4_twist_inverse(x,n,m,right_lo+m);
    //DO THE SAME FOR OTHER SIDE
    lut_tangent_8_inverse(x,m,right_lo);
    lut_cost_4_untwist_inverse(x,n,m,right_lo);

    //STITCH BACK TOGETHER (x^2n-i) and (x^2n+i) TO (x^4n+1)
    for (int i = right_lo; i < right_lo+m; ++i)
    {
      // printf("i = %d, m = %d\n",i,m );
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = -I*(temp - x[i+m]);
    }

    //Go from (x^4n-1) to (x^2n-1) and (x^2n+1)
    for(int i=lo; i < lo+m;++i){
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }
    m = m*2;

    //Go from (x^8n +1 to x^4n -1 and x^4n +1)
    for(int i=lo; i < lo+m;++i){
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }

  }
}

void lut_remove_base(double complex *x,int lo,int m, int n)
{
  int j =1;
  for (int i = lo+1; i < lo+m; ++i)
  {
    x[i] = x[i] / lut_s(n,j);
    ++j;
  }
}

void lut_tangent_4(double complex *x,int n, int lo)
{
    double complex temp;
    if(n == 2){
      // printf("LEVEL2\n");
      temp = x[lo];
      x[lo] = temp + x[lo+1];
      x[lo+1] = temp - x[lo+1];
    // print_complex(x,8);
    }
  else if(n>2)
  {
     int m = n/2;
      //Go from (x^4n +1 to x^2n -1 and x^2n +1)
      for(int i=lo; i < lo+m;++i){
        temp = x[i];
        x[i] = temp + x[i+m];
        x[i+m] = temp - x[i+m];
      }
      //Do recursive step for (x^2n -1)
      lut_tangent_4(x,m,lo);

      lo = lo+m;
      m = m/2;
      //Go from (x^2n +1 to x^n -i and x^n +i)
      for (int i = lo; i < lo+m; ++i)
      {
        // printf("i = %d, m = %d\n",i,m );
        temp = x[i];
        x[i] = temp + I * x[i+m];
        x[i+m] = temp - I * x[i+m];
      }
      lut_twist(x,n,m,lo);
      lut_base_change(x,m,-1,m,lo);
      lut_tangent_8(x,m,lo);
      lut_untwist(x,n,m,lo+m);
      lut_base_change(x,m,-1,m,lo+m);
      // printf("AFTER UNTWISTING\n");
      // print_complex(x,CPLXDIM);
      lut_tangent_8(x,m,lo+m);
  } 
}

void lut_tangent_4_inverse(double complex *x,int n, int lo)
{
  double complex temp;
    if(n == 2){
      // printf("LEVEL2\n");
      temp = x[lo];
      x[lo] = temp + x[lo+1];
      x[lo+1] = temp - x[lo+1];
    // print_complex(x,8);
    }
  else if(n > 2)
  {
      // printf("n = %d lo = %d\n",n,lo );
      int m = n/4;
      lo = lo+n/2;
      // printf("m = %d lo = %d\n",m,lo );
      lut_tangent_8_inverse(x,m,lo+m);
      lut_remove_base(x,lo+m,m,m);
      lut_twist(x,n,m,lo+m);
      
      lut_tangent_8_inverse(x,m,lo);
      lut_remove_base(x,lo,m,m);
      lut_untwist(x,n,m,lo);

      for (int i = lo; i < lo+m; ++i)
      {
        // printf("i = %d, m = %d\n",i,m );
        temp = x[i];
        x[i] = temp + x[i+m];
        x[i+m] = (temp - x[i+m])*-I;
      }
      // printf("AFTER\n");
      // print_complex(x,CPLXDIM);
      m = m*2;
      lo = lo -m;
      // printf("m = %d lo = %d\n",m,lo );
      lut_tangent_4_inverse(x,m,lo);
      for(int i=lo; i < lo+m;++i){
        temp = x[i];
        x[i] = temp + x[i+m];
        x[i+m] = temp - x[i+m];
      }
  }
}

void lut_tangent_forward(double complex *x)
{
  lut_twist(x,ROOTDIM,CPLXDIM,0);
  lut_tangent_4(x,CPLXDIM,0);
}

void lut_tangent_backward(double complex *x)
{
  lut_tangent_4_inverse(x,CPLXDIM,0);
  lut_untwist(x,ROOTDIM,CPLXDIM,0);
}