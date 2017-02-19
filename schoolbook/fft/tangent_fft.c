#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "tangent_fft.h"
#include "../support.h"
#include "../../test.h"

/******************************************************************
*
* Tangent FFT supporting functions
*
******************************************************************/
double s(int n, int k) 
{
  if (n <= 4) 
    return 1.0;

  int k4 = k % (n/4);

  if (k4 <= n/8) 
    return (s(n/4,k4) * cos(2.0 * M_PI * (double)k4 / (double)n));
  
  return (s(n/4,k4) * sin(2.0f * M_PI * (double)k4 / (double)n));
}

void base_change(double complex *x,int n,int n_new,int m,int lo)
{ int j = 1;
  if(n_new == -1)
  {
    for (int i = lo+1; i < lo+m; ++i)
    {
      x[i] = x[i] * s(n,j);
      ++j;
    }
  }
  else
  {
    double temp;
    for (int i = lo+1; i < lo+m; ++i)
    { 
      temp = s(n_new,j)/s(n,j);
      // printf("temp = %f, x[%d] = %f + I %f\n",temp,i,creal(x[i]),cimag(x[i]));
      x[i] = x[i] * temp;
      ++j;
    }
  }
}

void cost_4_twist(double complex *cplx_x,int n,int m,int lo)
{
  int j = 1;
  double complex temp;
  // printf("lo = %d m = %d n = %d\n",lo,m,n);
  for (int i = lo+1; i < lo+m; ++i)
  {
    temp = s(m,j)/s(n,j);
    // printf("temp = %f + i %f\n",creal(temp),cimag(temp));
    // printf("W(%d,%d) = %f + i %f\n",n,j,creal(W(n,j)),cimag(W(n,j)));
    // printf("temp = %f + i %f, x[%d] = %f + I %f\n",creal(temp),cimag(temp),i,creal(cplx_x[i]),cimag(cplx_x[i]));
    cplx_x[i] = cplx_x[i] * W(n,j) *temp;
    ++j;
  }
}

void cost_4_untwist(double complex *cplx_x,int n,int m,int lo)
{
  int j = 1;
  double complex temp;
  for (int i = lo+1; i < lo+m; ++i)
  {
    temp = s(m,j)/s(n,j);
    cplx_x[i] = cplx_x[i] * conj(W(n,j))*temp;
    ++j;
  }
}

void cost_4_untwist_inverse(double complex *cplx_x,int n,int m,int lo)
{
  int j = 1;
  double complex temp;
  for (int i = lo+1; i < lo+m; ++i)
  {
    temp = s(n,j)/s(m,j);
    cplx_x[i] = cplx_x[i] * conj(W(n,j))*temp;
    ++j;
  }
}

void remove_base(double complex *x,int lo,int m, int n)
{
  int j =1;
  for (int i = lo+1; i < lo+m; ++i)
  {
    x[i] = x[i] / s(n,j);
    ++j;
  }
}

void cost_4_twist_inverse(double complex *cplx_x,int n,int m,int lo)
{
  int j = 1;
  double complex temp;
  // printf("lo = %d m = %d n = %d\n",lo,m,n);
  for (int i = lo+1; i < lo+m; ++i)
  {
    temp = s(n,j)/s(m,j);
    // printf("temp = %f + i %f\n",creal(temp),cimag(temp));
    // printf("W(%d,%d) = %f + i %f\n",n,j,creal(W(n,j)),cimag(W(n,j)));
    // printf("temp = %f + i %f, x[%d] = %f + I %f\n",creal(temp),cimag(temp),i,creal(cplx_x[i]),cimag(cplx_x[i]));
    cplx_x[i] = cplx_x[i] * W(n,j) *temp;
    ++j;
  }
}

void tangent_8(double complex *x,int n,int lo)
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
    tangent_8(x,m,lo);
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
      cost_4_twist(x,n,m,right_lo);
      tangent_8(x,m,right_lo);
      //DO COST-4 UNTWIST
      cost_4_untwist(x,n,m,right_lo+m);
      tangent_8(x,m,right_lo+m);
      //BACK TO THE LEFT
      //NOW TO A BASE TWIST GO FROM X^2n-1 BASE S(8n,k) to X^2n-1 BASE S(2n,k)
      base_change(x,n,m,m,lo);
      tangent_8(x,m,lo);
      lo = lo+m;
      //FINISH THE LEFT_RIGHT BRANCH
      //ALSO BASE CHANGE FOR MID GO FROM X^2n+1 BASE S(8n,k) to X^2n+1 BASE S(4n,k)
      // printf("BASECHANGE!!!\n");
      // print_complex(x,8);
      base_change(x,n,n/2,m,lo);
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
      // base_change(x,lo,m,n/2,n/8);
      // twist(x,n/2,m,lo);
      cost_4_twist(x,n/2,m,lo);
      tangent_8(x,m,lo);
      // base_change(x,lo+m,m,n/2,n/8);
      // untwist(x,n/2,m,lo+m);
      cost_4_untwist(x,n/2,m,lo+m);
      tangent_8(x,m,lo+m);
  }
}

void tangent_8_inverse(double complex *x,int n,int lo)
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
    tangent_8(x,m,lo);
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
    tangent_8_inverse(x,m,lo+m);
    cost_4_twist_inverse(x,n/2,m,lo+m);
    tangent_8_inverse(x,m,lo);
    cost_4_untwist_inverse(x,n/2,m,lo);
    //GO FROM (X^n-i) AND (x^n+i) TO (X^2n+1)
    for (int i = lo; i < lo+m; ++i)
    {
      // printf("i = %d, m = %d\n",i,m );
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = -I*(temp - x[i+m]);
    }
    m = m*2;
    //NOW WE NEED TO CHANGE FROM BASE S(4n,k) TO s(8n,k)
    // printf("BASECHANGE!!!\n");
    // print_complex(x,8);
    base_change(x,n/2,n,m,lo);
    // printf("BASECHANGE!!!\n");
    // print_complex(x,8);
    lo = lo-m;
    //BACK TO THE LEFT WE FIRST NEED TO GET LEFT BRANCH (x^2n-1) BASE S(2n,k)
    tangent_8_inverse(x,m,lo);
    //NOW TO A BASE CHANGE GO FROM X^2n-1 BASE S(2n,k) to X^2n-1 BASE S(8n,k)
    base_change(x,n/4,n,m,lo);
    //NOW GO BACK TO THE RIGHT PART
    //GET RIGHT (x^2n-1)
    tangent_8_inverse(x,m,right_lo+m);
    //TWIST (x^2n-1) to (x^2n+i) with base S(8n,k)
    cost_4_twist_inverse(x,n,m,right_lo+m);
    //DO THE SAME FOR OTHER SIDE
    tangent_8_inverse(x,m,right_lo);
    cost_4_untwist_inverse(x,n,m,right_lo);

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

void tangent_4(double complex *x,int n, int lo)
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
      tangent_4(x,m,lo);

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
      twist(x,n,m,lo);
      base_change(x,m,-1,m,lo);
      tangent_8(x,m,lo);
      untwist(x,n,m,lo+m);
      base_change(x,m,-1,m,lo+m);
      // printf("AFTER UNTWISTING\n");
      // print_complex(x,CPLXDIM);
      tangent_8(x,m,lo+m);
  } 
}

void tangent_4_inverse(double complex *x,int n, int lo)
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
      tangent_8_inverse(x,m,lo+m);
      remove_base(x,lo+m,m,m);
      twist(x,n,m,lo+m);
      
      tangent_8_inverse(x,m,lo);
      remove_base(x,lo,m,m);
      untwist(x,n,m,lo);

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
      tangent_4_inverse(x,m,lo);
      for(int i=lo; i < lo+m;++i){
        temp = x[i];
        x[i] = temp + x[i+m];
        x[i+m] = temp - x[i+m];
      }
  }
}

void tangent_forward(double complex *x,int dim)
{
  tangent_4(x,dim,0);
}

void tangent_backward(double complex *x,int dim)
{
  tangent_4_inverse(x,dim,0);
}