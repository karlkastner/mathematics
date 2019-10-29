// Wed Oct 12 01:24:07 MSK 2011
// Karl Kästner, Berlin 

#ifndef SUMMATION_H
#define SUMMATION_H 1

#include <cstdint>
#include <cstring>
#include <ctime>
extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

#define half __float16
#define quad __float128

// TODO algorithm with sorting first

// trivial summation err ~ n*eps
template <class T>
T sum_serial( T const * const array, const unsigned int size )
{
	// trivial
	T sum = 0;
	for (unsigned int i=0; i<size; i++)
	{
		sum += array[i];
	}
	return sum;
} // sum_serial

// Kahan summation
// err ~ eps
template <class T>
T sum_compensated( T const * const array, const unsigned int size )
{
	T sum = 0.0;
	T remainder = 0.0;
	for (uint32_t i=0; i<size; i++)
	{
		T value = array[i] - remainder;
		T tmp = sum + value;
		remainder = (tmp - sum) - value;
		sum = tmp;
	}
	return sum;
} // sum_compensated

// element wise copy-cast
// B = (T2) A
template <class T1, class T2>
void copy_cast(T1 const * const A, T2 * const B, const uint32_t n)
{
	for (unsigned int i=0; i<n; i++)
	{
		B[i] = (T2) A[i];
	}
} // copy_cast

// elementwise product (see also: inner product, outer product)
// C = A.*B
template <class T>
void prod_ew(T const * const A, T const * const B, T * const C, const uint32_t n )
{
	for (unsigned int i=0; i<n; i++)
	{
		C[i] = A[i]*B[i];
	}
} // prod_ew

// pairwise accumulation
// err ~ log(n)*eps
// TODO implement cache oblivious version with reduced memory footprint in 32k chunks
template <class T>
T sum_pairwise(T const * const array, const int size )
{
	// save claculation, copy first
	T * const buffer = new T[size];
	memcpy(buffer, array, size*sizeof(T));
	/*int n = 1;
	while (n<size)
	{
		int m = 2*n;
		for(int i=0; i<size-n; i+=m)
		{
			buffer[i] += buffer[i+n];
		}
		n = m;
	}*/
// allocate buffer
/*
	sizte_t m = 32k/sizeof(T);
m > n // ?
	T* buffer = new T[n/m];
	for (int i=0; i<n/m; i++)
	{
		// recursive call
		buffer[i] = sum_P(&array, m);
	}
	todo last element
/ *
	for (int i=0; i<n_chunk-1; i++)
	{
		// copy chunk into auxillary array
		memcpy(chunk, array, CACHESIZE);
		// reduce chunk
		buffer[i] = reduce(chunk);
	}
	//
	last chunk
*/

	unsigned int n=size;
	while(n>1)
	{
		// add odd element
		if (0 != (n&1))
		{
			buffer[0] += buffer[n-1];
		}
		n = (n>>1);
		// accumulate regular elements
		for (unsigned int i=0; i<n; i++)
		{
			buffer[i] += buffer[i+n];	
		}
	}
	T sum = buffer[0];
	delete [] buffer;
	return sum;
} // sum_pairwise

template <class T>
void randn(T * const A, uint32_t n)
{
  	const gsl_rng_type * RT;
	gsl_rng * r;
	gsl_rng_env_setup();
  	RT = gsl_rng_default;
	r = gsl_rng_alloc (RT); //gsl_rng_mt19937)
	struct timespec tp;
	clock_gettime(CLOCK_REALTIME, &tp);
	gsl_rng_set (r, tp.tv_nsec);
	//gsl_rng_set (r, clock());//time(NULL));

	// todo seed

	for (uint32_t i=0; i<n; i++)
	{
		A[i] = (T) gsl_ran_ugaussian(r);
//		A[i] = (T) gsl_ran_gaussian(r, 1.0);
	}

	gsl_rng_free (r);
} // randn

#endif // ifndef SUMMATION_H

