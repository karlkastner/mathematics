// Karl Kästner, Berlin
// Wed Nov  9 01:28:45 MSK 2011
#include <cstdio>
#include <cstdlib>

extern "C" {
#include "f2c.h"
#include "clapack.h"
#include "cblas.h"
}

// golub 489
// Laplacian, once larger cpp: 2^20 => round off error or bug ?
// programm by yourself !!!
// SM,LA,SA 2^9-2^10, SM many more > 2^20
// matlab LA: 554, SM not reached
// m=8; n=2^25; D = sparse(-2*ones(n,1)); E = sparse(ones(n,1)); T=diag(D) + diag(E(1:end-1),-1) + diag(E(1:end-1),+1); tic, eigs(T,m,'SM'), toc
// D=full(D); E=full(E); save('e.dat','-ascii','E'); save('d.dat','-ascii','D')
// g++ eigs.cpp -lblas -L/usr/lib/atlas-base -llapack -Wall


extern "C" {
int dhseqr_( char * JOB, char * COMPZ, int * N, int * ILO, int * IHI, double * H,
            int * LDH, double * WR, double * WI, double * Z, int * LDZ,
            double * WORK, int * LWORK, int * INFO );
}

extern "C" {
int dstebz_(char *range, char *order, integer *n, double *vl,
double *vu, integer *il, integer *iu, double *abstol, double *d__, double *e,
integer *m, integer *nsplit, double *w, integer *iblock, integer *isplit, double
*work, integer *iwork, integer *info);
}

void error(char const * const name_cA, const int line_i)
{
	fprintf(stderr,"Error: %s %d\n",name_cA, line_i);
	exit(EXIT_FAILURE);
}

// read a vector from a file
int read(const char * const name_cA, double * const v_dA, const int n_i)
{
	FILE * fh = fopen(name_cA, "r");
	if (NULL == fh) error(__FILE__, __LINE__);
	for (int i=0; i<n_i; i++)
	{
		int ret = fscanf(fh,"%lf", &v_dA[i]);
		//printf("%f %d\n", v_dA[i],ret);
	}
	fclose(fh);
	return 0;
} // read

int main(int argc, char * argv[])
{
	if (argc < 5)
	{
		error(__FILE__, __LINE__);
	}

	// matrix rank
	int n_mat_i = atoi(argv[1]);
	// number of requested eigenvalues
	// also return value, number of eigenvalues found
	int n_eig_i = atoi(argv[2]);
	// diagonal data
	double * d_dA = new double[n_mat_i];			// 1 n
	read(argv[3], d_dA, n_mat_i);
	// sub-and-super-diagonal data
	double * s_dA = new double[n_mat_i];			// 2 n
	read(argv[4], s_dA, n_mat_i);
	// select eigenvalues by ordering
	char range_c = 'I';
	// sort eigenvalues
	char order_c = 'E';
	// upper and lower bounds (unused)
	double vl_d = NULL;
	double vu_d = NULL;
	// first and last eigenvalue
	int nl_i = n_mat_i - n_eig_i + 1;
	int nu_i = n_mat_i;
	double abstol_d = 1e-15;
	// number of diagonal blocks
	int n_split_i = 1;
	// eigenvalues, at least size M
	double * eig_dA = new double[n_mat_i];			// 3 n
	// return value, size N
	int * block_iA = new int[n_mat_i];			// 4 n
	// return value, size N
	int * split_iA = new int[n_mat_i];			// 5 n
	// workspace, size 4 n
	double * work_dA = new double[4*n_mat_i];		// 9 n
	// workspace, size 3 n
	int * work_iA = new int[3*n_mat_i];			// 12 n
	// return value
	int info_i;

	// find eigenvalues
	dstebz_( &range_c, &order_c, &n_mat_i, &vl_d, &vu_d, &nl_i, &nu_i, &abstol_d, \
		d_dA, s_dA, &n_eig_i, &n_split_i, eig_dA, block_iA, split_iA, \
			work_dA, work_iA, &info_i );
	// check return value
	if (0 != info_i)
	{
		error(__FILE__,__LINE__);
	}

	// print eigenvalues
	for (int i=0; i<n_eig_i; i++)
	{
		printf("%1.16e\n", eig_dA[i]);
	}

	// tidy up
	delete [] work_iA;
	delete [] work_dA;
	delete [] split_iA;
	delete [] block_iA;
	delete [] eig_dA;
	
	return EXIT_SUCCESS;
} // function main

