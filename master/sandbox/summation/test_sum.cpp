#include "summation.h"

#include <cstdlib>
#include <cmath>
#include <cstdio>

int main()
{
	int k=7;
	double sum;
	double err_serial;
	double err_pairwise;
	double err_compensated;

	for (int i=1; i<=k; i++)
	{
		uint32_t n = (uint32_t) pow(10.0,i);
		// allocate space
		double *A = new double[n];
		double *B = new double[n];
		double *C = new double[n];
		float  *C_single = new float[n];
		
		// generate two random input arrays
		randn(A, n);
		randn(B, n);

		// inner product
		prod_ew(A, B, C, n);
		// copy to single precission
		copy_cast(C, C_single, n);
		sum	        = sum_serial(C, n);
		err_serial      = fabs(sum_serial(C_single, n) - sum);
		err_pairwise    = fabs(sum_pairwise(C_single,n) - sum);
		err_compensated = fabs(sum_compensated(C_single,n) - sum);

		printf("%8u %e %e %e %e ", n, sum, err_serial, err_pairwise, err_compensated);

		// norm
		prod_ew(A, A, C, n);
		// copy to single precission
		copy_cast(C, C_single, n);
		
		sum 	        = sum_serial(C, n);
		err_serial      = fabs(sum_serial(C_single, n) - sum);
		err_pairwise    = fabs(sum_pairwise(C_single,n) - sum);
		err_compensated = fabs(sum_compensated(C_single,n) - sum);
		
		printf("%e %e %e %e\n", sum, err_serial, err_pairwise, err_compensated);

		delete [] A;
		delete [] B;
		delete [] C;
		delete [] C_single;
	}

	return EXIT_SUCCESS;
} // main

