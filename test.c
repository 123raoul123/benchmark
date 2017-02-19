#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "test.h"
#include "schoolbook/zeropad_mul.h"
#include "schoolbook/twisting_mul.h"
#include "lut/lut_mul.h"
#include "vector/vec_mul.h"

#define NTESTS 1000
#define NALGO 21
#define THRESHOLD 2
#define DEBUG 0


uint32_t myabs(uint32_t a)
{
  return (a>(1UL<<31)?-a:a);
}

int compare(const void * elem1, const void * elem2)
{	
	uint64_t x = *((uint64_t*)elem1);
	uint64_t y = *((uint64_t*)elem2);
	if(x > y) 
		return 1;
	else if(x < y) 
		return -1;
	return 0;
}

uint64_t rdtsc()
{
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

void cycle_meassure(){
  FILE *urandom = fopen("/dev/urandom", "r");
  ring_t re,r,x,y;
  int n,i;
  uint64_t start,end; 
  uint64_t cycles[NTESTS];
  
  lut_init();
  vec_init();

  for (int j = 0; j < NALGO; ++j)
  {
  	int success = 0;
    for(n=0;n<NTESTS;n++)
    { 

		for(i=0;i<REALDIM;++i){
		x.v[i] = fgetc(urandom);
		y.v[i] = fgetc(urandom);
		}
		switch(j)
		{
			case 0:
			  	start = rdtsc();
			    schoolbook(&r,&x,&y);
			    end = rdtsc();
			    break;
			case 1:
			    start = rdtsc();
			    zeropad_fiduccia(&r,&x,&y);
			    end = rdtsc();
			    break;
			case 2:
			    start = rdtsc();
				zeropad_twisted(&r, &x, &y);
				end = rdtsc();
				break;
			case 3:
			    start = rdtsc();
			    zeropad_split_radix(&r,&x,&y);
			    end = rdtsc();
			    break;
			case 4:
			    start = rdtsc();
				zeropad_tangent(&r,&x,&y);
				end = rdtsc();
				break;
			case 5:
				start = rdtsc();
				twist_fiduccia(&r,&x,&y);
				end = rdtsc();
				break;
			case 6:
				start = rdtsc();
				twist_negacyclic(&r,&x,&y);
				end = rdtsc();
				break;
			case 7:
				start = rdtsc();
				twist_twisted(&r,&x,&y);
				end = rdtsc();
				break;
			case 8:
				start = rdtsc();
				twist_split_radix(&r,&x,&y);
				end = rdtsc();
				break;
			case 9:
				start = rdtsc();
				twist_tangent(&r,&x,&y);
				end = rdtsc();
				break;
			case 10:
				start = rdtsc();
				lut_negacyclic(&r,&x,&y);
				end = rdtsc();
				break;
			case 11:
				start = rdtsc();
				lut_split_radix(&r,&x,&y);
				end = rdtsc();
				break;
			case 12:
				start = rdtsc();
				lut_tangent(&r,&x,&y);
				end = rdtsc();
				break;
			case 13:
				start = rdtsc();
			    simd_negacyclic(&r,&x,&y);
			    end = rdtsc();
			case 14:
			  	start = rdtsc();
			    simd_nonrec_negacyclic(&r,&x,&y);
			    end = rdtsc();
			case 15:
			  	start = rdtsc();
			    two_layer_negacyclic(&r,&x,&y);
			    end = rdtsc();
			    break;
			case 16:
			  	start = rdtsc();
			    three_layer_negacyclic(&r,&x,&y);
			    end = rdtsc();
			    break;
			case 17:
			    start = rdtsc();
			    sr_vector_mul(&r,&x,&y);
			    end = rdtsc();
			    break;
			case 18:
			    start = rdtsc();
				sr_vector_nonrec_mul(&r, &x, &y);
				end = rdtsc();
				break;
			case 19:
			    start = rdtsc();
			    fftw_mul(&r,&x,&y);
			    end = rdtsc();
			    break;
			case 20:
			    start = rdtsc();
				fftw_nega_mul(&r,&x,&y);
				end = rdtsc();
				break;

		}
		if(DEBUG)
		{
			schoolbook(&re,&x,&y);
			bool error = false;
	    	for(i=0;i<REALDIM;i++)
		    {

		      if(myabs(r.v[i] - re.v[i]) > THRESHOLD){
		        printf("school: %u \n", re.v[i]);
		        printf("test: %u \n",r.v[i]);
		        printf("difference: %u\n\n",myabs(r.v[i] - re.v[i]));
		        error = true;
		      }
		      
		    }
		    if(!error)
		      success++;
		}
		cycles[n] = end - start;
    }
    if(DEBUG)
    	printf("Ammount of successful multiplications: %d\n", success);

    qsort(cycles,sizeof(cycles)/sizeof(*cycles),sizeof(*cycles),compare);
    switch(j)
    {
    	case 0:
    		printf("ZEROPADDING\n");
    		printf("Schoolbook & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
    		break;
		case 1:
			printf("Normal FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
		case 2:
			printf("Twisted FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
		case 3:
			printf("Split-Radix FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
		case 4:
			printf("Tangent FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
    	case 5:
    		printf("TWIST\n");
    		printf("FFT& %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
    		break;
		case 6:
			printf("Negacyclic FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
		case 7:
			printf("Twisted FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
		case 8:
			printf("Split-Radix FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
		case 9:
			printf("Tangent FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
		case 10:
			printf("LOOKUPTABLES\n");
			printf("Negacyclic FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
		case 11:
			printf("Split-Radix FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
		case 12:
			printf("Tangent FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
		case 13:
	    	printf("VECTORIZED\n");
      		printf("simd recursive Negacyclic FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
      		break;
	    case 14:
	    	printf("simd nonrecursive Negacyclic FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
      		break;
	    case 15:
	    	printf("TWO LAYER Negacyclic FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
      		break;
	    case 16:
	    	printf("THREE LAYER Negacyclic FFT & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
      		break;
	    case 17:
			printf("SplitRadix & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
	    case 18:
			printf("Non recursive SR & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
	    case 19:
			printf("FFTW & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;
	    case 20:
			printf("Twisted FFTW & %llu & %llu & %llu\\\\\n",cycles[249],cycles[499],cycles[749]);
			break;

    }
      
  }

  fclose(urandom); 
}

int main()
{
  cycle_meassure();

  return 0;
}