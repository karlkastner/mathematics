// Thu Jun  7 20:25:01 MSK 2012
// Karl Kästner, Berlin
import java.lang.Math.*;
import Jama.*;

public class Assemble_2d_phi_phi
{
	int MAX_NUM_THREADS = 8;
	int n_thread;
	Worker[] t;
	int n_first[];
	int n_last[];
	
	double[][] P;
	double [][][] P_local;
	int[][] T;
	double [][][] Phi;
	Potential_2D func;
	double [] w;
	double [][] b;
	int flag;
	int mb0;
	int nb0;
	int nt1;
	int nt2;
	int nv;
	double [][] buf;
	double [] determinant;

	// default constructor
	public Assemble_2d_phi_phi()
	{
		t       = new Worker[MAX_NUM_THREADS];
		n_first = new int[MAX_NUM_THREADS];
		n_last  = new int[MAX_NUM_THREADS];
	};

	public double [][] assemble(Mesh_2d mesh, Potential_2D func, double [] w, double [][] b, int flag) throws java.lang.InterruptedException
	{
		// accept arguments
		this.P    = mesh.P;
		this.P_local = mesh.P_local;
		this.T    = mesh.T;
		this.Phi  = mesh.Phi;
		this.determinant = mesh.determinant;
		this.func = func;
		this.w    = w;
		this.b    = b;
		this.flag = flag;
		nt1 = T.length;
		nt2 = T[0].length; // TODO dangerous, if T.length == 0

		// rounding not required in exact arithmetic
		nv = (int) Math.round(-1.5 + java.lang.Math.sqrt(2*nt2 + 0.25));

		// allocate memory
		if (0 != flag && null == func)
		{
			//  diagonal lumped mass matrix
			// test functions phi_i are chosen that integral is 1 at the points p_i and zero at the other points
			mb0 = nt2;
		} else {
			// non diagonal mass matrix
			mb0 = nt2*(nt2+1)/2;
		}
		nb0 = nt1*mb0;

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
	} // assemble

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

	// TODO instead of using two dimensional arrays use one dimensional arrays with stride (column wise like in FORTRAN)
	private void assemble_(int i)
	{
		// preallocate thread local buffers
		double[][] A    = new double[3][3];
		double[]   wfa  = new double[w.length];
		double[][] Va   = new double[nt2][(nv+1)*(nv+2)/2];
		double[][] Vq   = new double[w.length][(nv+1)*(nv+2)/2];
		Matrix mA       = new Matrix(A);
		Matrix mB       = new Matrix(b);
		Matrix mVq      = new Matrix(Vq);

		Matrix mVa = new Matrix(Va);
		double[][] A_   = new double[nt2][2];
		Matrix mC;

		double area;
		// integrate over each triangle
		for (int idx=n_first[i]; idx<=n_last[i]; idx++)
		{
			// fetch corner points
			if (null == P_local)
			{
				// no values are prefetched
				A[0][0] = 1; A[0][1] = P[T[idx][0]-1][0]; A[0][2] = P[T[idx][0]-1][1];
		                A[1][0] = 1; A[1][1] = P[T[idx][1]-1][0]; A[1][2] = P[T[idx][1]-1][1];
		                A[2][0] = 1; A[2][1] = P[T[idx][2]-1][0]; A[2][2] = P[T[idx][2]-1][1];
			} else {
				// values are prefetched
				A[0][0] = 1; A[0][1] = P_local[idx][0][0]; A[0][2] = P_local[idx][0][1];
				A[1][0] = 1; A[1][1] = P_local[idx][1][0]; A[1][2] = P_local[idx][1][1];
				A[2][0] = 1; A[2][1] = P_local[idx][2][0]; A[2][2] = P_local[idx][2][1];		    
			}

			// determine the triangle area
			if (null == determinant)
			{
				area = 0.5*Math.abs(mA.det());
			} else {
				area = 0.5*Math.abs(determinant[idx]);
			}
	
			// quadrature points in homogeneous Cartesian coordinates
			Matrix mQ = mB.times(mA);
	
			// no prefetch
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
				// 0, 2 changed to 0, lt2-1
				FEM.vander_2d(Va, A_, nv);
		
				// compute the coefficients of the polynomial test functions
				// C(i,:) * A(:,i) = 1; C(i,:) * A(:,j) = 0, i<>j
				mC = mVa.inverse();
			} else {
				// load precomputed test function coefficients
				mC = new Matrix(Phi[idx]);
			}
	
			// evaluate test function at the integration points
			// TODO these values should be constant and can be precomputed
			FEM.vander_2d(Vq, mQ.getMatrix(0, w.length-1, 1, 2).getArray(), nv);
	
			// phi is not perfectly symmetric!!!
			Matrix mPhi = mVq.times(mC);
			double[][] phi = mPhi.getArray();
	
			// evaluate the coefficient function at integration points
			if (null != func)
			{
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
	
			// stiffness matrix contribution
			// for all testfunction being 1 at point adx
			int m = mb0*idx;
			for (int adx=1; adx <= nt2; adx++)
			{
				// diagonal entry integral approximation
				double I = 0;
				for (int k=0; k<w.length; k++)
				{
					I = I + wfa[k]*phi[k][adx-1]*phi[k][adx-1];
				}
				buf[m][0] = T[idx][adx-1];
				buf[m][1] = T[idx][adx-1];
				buf[m][2] = 0.5*I;
				m = m+1;
				// off diagonal entries
				// exploit symmetry A(i,j) = A(j,i)
				if (0 == flag || null != func)
				{
					for (int bdx=(adx+1); bdx<=nt2; bdx++)
					{
						// integral approximation
						I = 0;
						for (int k=0; k<w.length; k++)
						{
							I = I + wfa[k]*phi[k][adx-1]*phi[k][bdx-1];
						}
						buf[m][0] = T[idx][adx-1];
						buf[m][1] = T[idx][bdx-1];
						buf[m][2] = I;
						m = m+1;
					} // for bdx
				}
			} // for adx
		} // for idx
	} // assemble_
} // class Assemble_2d_phi_phi

