CC=/usr/bin/gcc
CFLAGS=-Wall -Wextra -std=c11 -O3 -mavx2 -march=native
#DEBUGFLAGS = -Wall -Wextra -std=c11 -mavx2 -march=native -v -da -Q -g -O0
LIBS=-lm -L/usr/local/lib -lfftw3
OBJS_ZERO = zeropad_mul.o fiduccia.o split_radix_fft.o twisted_fft.o tangent_fft.o support.o twisting_mul.o
OBJS_LUT = lut_mul.o lut_negacyclic.o lut_split_radix.o lut_tangent.o
OBJS_VEC = vec_mul.o 2_layer_negacyclic.o sr_vector.o fftw.o sr_vec_nonrec.o simd_negacyclic.o 3_layer_negacyclic.o simd_nonrec_negacyclic.o simd_nonmerge.o nonrec_nonmerge.o
OBJS = $(OBJS_ZERO) $(OBJS_LUT) $(OBJS_VEC) glob_support.o

output: test.o
	$(CC) $(CFLAGS) $(LIBS) test.o $(OBJS) -o $@

test.o: test.c test.h $(OBJS)
	$(CC) $(CFLAGS) -c test.c

zeropad_mul.o: test.h schoolbook/zeropad_mul.c schoolbook/zeropad_mul.h fiduccia.o split_radix_fft.o twisted_fft.o tangent_fft.o
	$(CC) $(CFLAGS) -c schoolbook/zeropad_mul.c

twisting_mul.o: test.h schoolbook/twisting_mul.c schoolbook/twisting_mul.h fiduccia.o split_radix_fft.o twisted_fft.o tangent_fft.o support.o glob_support.o
	$(CC) $(CFLAGS) -c schoolbook/twisting_mul.c

vec_mul.o: test.h sr_vector.o fftw.o sr_vec_nonrec.o 2_layer_negacyclic.o simd_negacyclic.o 3_layer_negacyclic.o simd_nonrec_negacyclic.o simd_nonmerge.o nonrec_nonmerge.o
	$(CC) $(CFLAGS) -c vector/vec_mul.c

sr_vector.o: test.h vector/fft/sr_vector.c vector/fft/sr_vector.h
	$(CC) $(CFLAGS) -c vector/fft/sr_vector.c

2_layer_negacyclic.o: test.h vector/fft/2_layer_negacyclic.c vector/fft/2_layer_negacyclic.h
	$(CC) $(CFLAGS) -c vector/fft/2_layer_negacyclic.c

3_layer_negacyclic.o: test.h vector/fft/3_layer_negacyclic.c vector/fft/3_layer_negacyclic.h
	 $(CC) $(CFLAGS) -c vector/fft/3_layer_negacyclic.c

simd_negacyclic.o: test.h vector/fft/simd_negacyclic.c vector/fft/simd_negacyclic.h
	$(CC) $(CFLAGS) -c vector/fft/simd_negacyclic.c

simd_nonrec_negacyclic.o: test.h vector/fft/simd_nonrec_negacyclic.c vector/fft/simd_nonrec_negacyclic.h
	$(CC) $(CFLAGS) -c vector/fft/simd_nonrec_negacyclic.c

sr_vec_nonrec.o: test.h vector/fft/sr_vec_nonrec.c vector/fft/sr_vec_nonrec.h
	$(CC) $(CFLAGS) -c vector/fft/sr_vec_nonrec.c

simd_nonmerge.o: test.h vector/fft/simd_nonmerge.c vector/fft/simd_nonmerge.h
	$(CC) $(CFLAGS) -c vector/fft/simd_nonmerge.c

nonrec_nonmerge.o: test.h vector/fft/nonrec_nonmerge.c vector/fft/nonrec_nonmerge.h
	$(CC) $(CFLAGS) -c vector/fft/nonrec_nonmerge.c

fftw.o: test.h vector/fft/fftw.c vector/fft/fftw.h
	$(CC) $(CFLAGS) -I/usr/local/include -c vector/fft/fftw.c

lut_mul.o: test.h lut/lut_mul.h lut/lut_mul.c lut_negacyclic.o lut_tangent.o lut_split_radix.o
	$(CC) $(CFLAGS) -c lut/lut_mul.c

lut_negacyclic.o: test.h lut/fft/lut_negacyclic.c lut/fft/lut_negacyclic.h
	$(CC) $(CFLAGS) -c lut/fft/lut_negacyclic.c

lut_split_radix.o: test.h lut/fft/lut_split_radix.c lut/fft/lut_split_radix.h
	$(CC) $(CFLAGS) -c lut/fft/lut_split_radix.c

lut_tangent.o: test.h lut/fft/lut_tangent.c lut/fft/lut_tangent.h
	$(CC) $(CFLAGS) -c lut/fft/lut_tangent.c

glob_support.o: test.h glob_support.h glob_support.c
	$(CC) $(CFLAGS) -c glob_support.c

fiduccia.o: schoolbook/fft/fiduccia.h schoolbook/fft/fiduccia.c
	$(CC) $(CFLAGS) -c schoolbook/fft/fiduccia.c

split_radix_fft.o: schoolbook/fft/split_radix_fft.h schoolbook/fft/split_radix_fft.c support.o
	$(CC) $(CFLAGS) -c schoolbook/fft/split_radix_fft.c

twisted_fft.o: schoolbook/fft/twisted_fft.h schoolbook/fft/twisted_fft.c support.o
	$(CC) $(CFLAGS) -c schoolbook/fft/twisted_fft.c

tangent_fft.o: schoolbook/fft/tangent_fft.h schoolbook/fft/tangent_fft.c support.o
	$(CC) $(CFLAGS) -c schoolbook/fft/tangent_fft.c

support.o: schoolbook/support.h schoolbook/support.c test.h
	$(CC) $(CFLAGS) -c schoolbook/support.c

clean:
	-rm *.o output
