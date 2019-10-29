// Thu Jun  7 21:27:01 MSK 2012
// Karl Kästner, Berlin

import java.lang.Math.*;
import Jama.*;

public class Assemble_2d_dphi_dphi
{
	int MAX_NUM_THREADS = 8;
	int n_thread;
	Worker[] t; // = new Worker[MAX_NUM_THREADS];
	int n_first[]; // = new int[MAX_NUM_THREADS];
	int n_last[]; // = new int[MAX_NUM_THREADS];
	
	double[][] P;
	int[][] T;
	double [][][] P_local;
	double [][][] Phi;
	Potential_2D func;
	double [] w;
	double [][] b;
	double [] determinant;
	int flag;
	int mb0;
	int nb0;
	int nt1;
	int nt2;
	int nv;
	//double [][] buf;
	volatile double [][] buf;
	//final double [][] buf;

	public Assemble_2d_dphi_dphi()
	{
		t       = new Worker[MAX_NUM_THREADS];
		n_first = new int[MAX_NUM_THREADS];
		n_last  = new int[MAX_NUM_THREADS];
	};

	public double [][] assemble(Mesh_2d mesh, Potential_2D func, double[] w, double[][] b) throws java.lang.InterruptedException
	{
		// accept arguments
		this.P = mesh.P;
		this.T = mesh.T;
		this.Phi = mesh.Phi;
		this.determinant = mesh.determinant;
		this.func = func;
		this.w = w;
		this.b = b;
		this.P_local = mesh.P_local;
		nt1 = T.length;
		nt2 = T[0].length; // TODO dangerous, if T.length == 0

		// preallocate buffer
		nt1 = T.length;
		nt2 = T[0].length; //TODO dangerous if T is empty
		// rounding in exact arithmetic unnecessary
		nv = (int) Math.round(-1.5 + Math.sqrt(2*nt2 + 0.25));
		mb0 = nt2*(nt2+1)/2;
		nb0 = nt1*mb0;
	
		// allocate memory
		buf = new double[nb0][3];

		// query the number of requested threads
		String omp_n_thread_s = System.getenv("OMP_NUM_THREADS");
		int omp_n_thread;
		try {
			omp_n_thread = Integer.parseInt( omp_n_thread_s );
		} catch ( Exception e )
		{
			omp_n_thread = java.lang.Integer.MAX_VALUE;
		}

		// query the number of available CPUs
		int available_processors = Math.max(1, java.lang.Runtime.getRuntime().availableProcessors());

		// determine the number of threads
		n_thread = Math.min(Math.min( MAX_NUM_THREADS, available_processors), omp_n_thread);

		// Why do the workers have to be created twice?
		int n_slice = nt1/n_thread;
		// first thread variables
		t[0]       = new Worker();
		n_first[0] = 0;
		n_last[0]  = n_slice-1;
		for (int i=1; i<n_thread; i++)
		{	
			n_first[i] = n_first[i-1]+n_slice;
			n_last[i]  = n_last[i-1]+n_slice;
			t[i]       = new Worker();
		}
		// let the last thread process the remainder
		n_last[n_thread-1] = nt1-1;

		// start the threads
		for (int i=0; i<n_thread; i++)
		{
			t[i].start();
		}
	
		// join the threads
		for (int i=0; i<n_thread; i++)
		{
			t[i].join();
		}

		// return the results
		return buf;
	} // assemble()

	private class Worker extends Thread
	{
		public synchronized void run()
		{
			for (int i=0; i<n_thread; i++)
			{
				if (t[i] == this)
				{
					assemble_(i);
				}
			}
		} // run()
	} // class Worker

	private void assemble_(int i)
	{
	double[][] A    = new double[3][3];
	double[][] A_   = new double[nt2][2];
	double[][] Vq   = new double[w.length][(nv)*(nv+1)/2];
	double[][] dC_x = new double[(nv)*(nv+1)/2][(nv+1)*(nv+2)/2];
	double[][] dC_y = new double[(nv)*(nv+1)/2][(nv+1)*(nv+2)/2];
	double wfa[]    = new double[w.length];
	double [][] C;

	Matrix mA    = new Matrix(A);
	Matrix mB    = new Matrix(b);
	Matrix mVq   = new Matrix(Vq);
	Matrix mDC_x = new Matrix(dC_x);
	Matrix mDC_y = new Matrix(dC_y);

	double[][] Va   = new double[nt2][(nv+1)*(nv+2)/2];
	Matrix mVa   = new Matrix(Va);

	//System.out.println(i + " started" + n_first[i] + " " + n_last[i] + " " + nt1);

	// integrate over each triangle
	for (int idx=n_first[i]; idx<=n_last[i]; idx++)
	//for (int idx=0; idx<nt1; idx++)
	{
		// fetch corner points - get those points in prefetch routine
		if (null == P_local)
		{
			A[0][0] = 1; A[0][1] = P[T[idx][0]-1][0]; A[0][2] = P[T[idx][0]-1][1];
                	A[1][0] = 1; A[1][1] = P[T[idx][1]-1][0]; A[1][2] = P[T[idx][1]-1][1];
                	A[2][0] = 1; A[2][1] = P[T[idx][2]-1][0]; A[2][2] = P[T[idx][2]-1][1];
		} else {
			A[0][0] = 1; A[0][1] = P_local[idx][0][0]; A[0][2] = P_local[idx][0][1];
			A[1][0] = 1; A[1][1] = P_local[idx][1][0]; A[1][2] = P_local[idx][1][1];
			A[2][0] = 1; A[2][1] = P_local[idx][2][0]; A[2][2] = P_local[idx][2][1];
		}

		double area;

		// determine the triangle area
		if ( null == determinant)
		{
			area = 0.5*Math.abs(mA.det());
		} else {
			area = 0.5*Math.abs(determinant[idx]);
		}

		// integration points in homogeneous Cartesian coordinates
		Matrix mQ = mB.times(mA);

		if (null == Phi)
		{
			// fetch all point coordinates, not just corner point coordinates
			// TODO, these coordinates can be computed on the fly (except on curved boundaries)
			for (int j=0; j<nt2; j++)
			{
				A_[j][0] = P[T[idx][j]-1][0];
				A_[j][1] = P[T[idx][j]-1][1];
			}
			// construct the triangle point Vandermonde matrix spanning the test function polynomials
			// structure: 1 x y xy x^2 y^2 x^2y xy^2 x^3 y^3
			FEM.vander_2d(Va, A_, nv);
			//FEM.vander_2d(Va, mA.getMatrix(0, lt2-1, 1, 2).getArray(), nv);

			// compute test functions
			// C(i,:) * A(:,i) = 1; C(i,:) * A(:,j) = 0, i<>j
			Matrix mC = mVa.inverse();
			C = mC.getArray();
		} else {
			// load precomputed test functions
			C = Phi[idx];
		}

		// form first derivative of the test functions
		// structure: dphi([1 x y]') : [dc00 dc01 dc10]*[1 x y]'
		FEM.derivative_2d(dC_x, dC_y, C, nv);

		// evaluate test function derivative at the integration points
		// test function derivative evaluated at the points
		// rows: function, columns: points
		FEM.vander_2d(Vq, mQ.getMatrix(0, w.length-1, 1, 2).getArray(), nv-1);
		Matrix mDphi_x = mVq.times(mDC_x);
		Matrix mDphi_y = mVq.times(mDC_y);
		double [][] dphi_x = mDphi_x.getArray();
		double [][] dphi_y = mDphi_y.getArray();

		// evaluate the coefficient function at integration points
		if (null != func)
		{
			//f = feval(func,q(:,2:3));
			double[] f = func.func(mQ.getMatrix(0, w.length-1, 1, 2).getArray());

			// premultiply values
			for (int k=0; k<w.length; k++)
			{
				wfa[k] = w[k]*f[k]*area;
			}
		} else {
			// premultiply values
			for (int k=0; k<w.length; k++)
			{
				wfa[k] = w[k]*area;
			}
		}

		// mass matrix contributions
		// for all testfunction being 1 at point adx
		int m = mb0*idx;
		for (int adx=1; adx<=nt2; adx++)
		{
			// diagonal entry integral approximation
			double I = 0;
			for (int k=0; k<w.length; k++)
			{
				I = I + wfa[k]*( dphi_x[k][adx-1]*dphi_x[k][adx-1]
					+ dphi_y[k][adx-1]*dphi_y[k][adx-1] );
			}
			buf[m][0] = T[idx][adx-1];
			buf[m][1] = T[idx][adx-1];
			buf[m][2] = -0.5*I;
			m = m+1;
			// off diagonal entries
			// exploit symmetry A(i,j) = A(j,i)
			// TODO, not necessary if lumped diagonal matrix flag is set
			for (int bdx=(adx+1); bdx<=nt2; bdx++)
			{
				// integral approximation
				I = 0;
				for (int k=0; k<w.length; k++)
				{
					I = I + wfa[k]*( dphi_x[k][adx-1]*dphi_x[k][bdx-1]
						+ dphi_y[k][adx-1]*dphi_y[k][bdx-1] );
				}
				buf[m][0] = T[idx][adx-1];
				buf[m][1] = T[idx][bdx-1];
				buf[m][2] = -I;
				m = m+1;
			} // for bdx
		} // for adx
	} // for idx
	//return buf;
	//System.out.println(i + " finnished"); // + n_first[i] + " " + n_last[i] + " " + nt1);
} // assemble_()
} // class Assemble_2d_dphi_dphi


