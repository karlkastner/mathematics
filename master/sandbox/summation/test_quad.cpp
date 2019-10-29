// Sep 12  2010
// Karl Kästner, Berlin

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstring>
// todo timing
// Loss of significance - auslöschung

const int MAX = 10*1024*1024;


int main()
{
	time_t t;
	time(&t);
	srand(t);

	//half
	float  *f32_a  = (float*) calloc(MAX, sizeof(float));
	double *f64_a  = (double*) calloc(MAX, sizeof(double));
	quad   *f128_a = (quad*) calloc(MAX, sizeof(quad));
	
	// generate random numbers
	for (int i=0; i<MAX; i++)
	{
		// RAND_MAX is 31 bit on my system
		// so four rands fill roughly a quad
		// its normed to fit into double and float
		quad q = ( (quad) rand()*rand()*rand()*rand())/((quad) rand()*rand()*rand()*rand());
//		quad q = (quad) 1.0 / ( (quad) rand()*pow(2,96) + (quad) rand()*pow(2,64) + (quad) rand()*pow(2,32) + (quad) rand() );
		// quad q = (quad) ( (quad) rand()*pow(2,96) + (quad) rand()*pow(2,64) + (quad) rand()*pow(2,32) + (quad) rand() ) /
		// ( (quad) rand()*pow(2,96) + (quad) rand()*pow(2,64) + (quad) rand()*pow(2,32) + (quad) rand() );
		f128_a[i] = q;
		f64_a[i] = (double) q;
		f32_a[i] = (float) q;
	}
	

	// sum trivial
	float  sumT_32 = sumT(f32_a, MAX);
	double sumT_64 = sumT(f64_a, MAX);
	quad   sumT_128 = sumT(f128_a, MAX);

	// sum Kahan
	float  sumK_32 = sumK(f32_a, MAX);
	double sumK_64 = sumK(f64_a, MAX);
	quad   sumK_128 = sumK(f128_a, MAX);

	// sum pairwise
	float  sumP_32 = sumP(f32_a, MAX);
	double sumP_64 = sumP(f64_a, MAX);
	quad   sumP_128 = sumP(f128_a, MAX);

	double errK_32 = fabs( (double) ( sumK_128 - (quad) sumK_32));
	double errK_64 = fabs( (double) ( sumK_128 - (quad) sumK_64));
	double errK_128 = fabs( (double) ( sumK_128 - (quad) sumK_128)); // test

	double errP_32 = fabs( (double) ( sumK_128 - (quad) sumP_32));
	double errP_64 = fabs( (double) ( sumK_128 - (quad) sumP_64));
	double errP_128 = fabs( (double) ( sumK_128 - (quad) sumP_128));

	double errT_32 = fabs( (double) ( sumK_128 - (quad) sumT_32));
	double errT_64 = fabs( (double) ( sumK_128 - (quad) sumT_64));
	double errT_128 = fabs( (double) ( sumK_128 - (quad) sumT_128));

	double sum = (double)sumK_128;
	std::cout<<(double)sumK_128<<" "<<(double)sumP_128<<" "<<(double)sumT_128<<std::endl;

	std::cout<<std::scientific<<sizeof(quad)<<" "<<std::scientific<<sizeof(double)<<" "<<std::scientific<<sizeof(float)<<std::endl;
	std::cout<<std::scientific<<(double)errK_128<<" "<<std::scientific<<log((sum/errK_128))/log(2.0)<<std::endl;
	std::cout<<std::scientific<<(double)errP_128<<" "<<std::scientific<<log((sum/errP_128))/log(2.0)<<std::endl;
	std::cout<<std::scientific<<(double)errT_128<<" "<<std::scientific<<log((sum/errT_128))/log(2.0)<<std::endl;
	std::cout<<std::scientific<<(double)errK_64<<" "<<std::scientific<<log((sum/errK_64))/log(2.0)<<std::endl;
	std::cout<<std::scientific<<(double)errP_64<<" "<<std::scientific<<log((sum/errP_64))/log(2.0)<<std::endl;
	std::cout<<std::scientific<<(double)errT_64<<" "<<std::scientific<<log((sum/errT_64))/log(2.0)<<std::endl;
	std::cout<<std::scientific<<(double)errK_32<<" "<<std::scientific<<log((sum/errK_32))/log(2.0)<<std::endl;
	std::cout<<std::scientific<<(double)errP_32<<" "<<std::scientific<<log((sum/errP_32))/log(2.0)<<std::endl;
	std::cout<<std::scientific<<(double)errT_32<<" "<<std::scientific<<log((sum/errT_32))/log(2.0)<<std::endl;
}

