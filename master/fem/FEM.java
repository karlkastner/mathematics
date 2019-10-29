// Sat Jun  2 23:53:25 MSK 20k2
// Karl Kästner, Berlin
import java.lang.Math.*;
import java.util.*;
import Jama.*;

public final class FEM
{
	// create the Vandermonde matrix of powers in 2D
	// const * const and * const in java ???
	public static final double [][] vander_2d(double[][] V, final double[][] p, final int n)
	{
		// todo restructure, that it works without pp-buffers (first fill in 1,x,x^2,y,y^2, .. then mixed
		// precalculate powers of P here
		double[][] pp0 = new double[p.length][n+1];
		double[][] pp1 = new double[p.length][n+1];
		for (int k=0; k<p.length; k++)
		{
			pp0[k][0] = 1;
			pp1[k][0] = 1;
		}
		for (int i=1; i<=n; i++)
		{
			for (int k=0; k<p.length; k++)
			{
				pp0[k][i] = pp0[k][i-1]*p[k][0];
				pp1[k][i] = pp1[k][i-1]*p[k][1];
			}
		}
		// set up the vandermonde matrix
		int m = 0;
		for (int idx=0; idx<=n; idx++)
		{
			for (int jdx=0; jdx<=idx; jdx++)
			{
				for (int k=0; k<p.length; k++)
				{
					// just access powers of P here
					//V[k][m] = java.lang.Math.pow(p[k][0],idx-jdx) * java.lang.Math.pow(p[k][1],jdx);
					V[k][m] = pp0[k][idx-jdx] * pp1[k][jdx];
				} // for k
				m = m+1;
			} // for jdx
		} // for idx
		return V;
	} // vander_2d

	// Va has to be preallocated sufficiently
	public static final double [][] vander_3d(double [][] V, final double [][] X, final int p)
	{
		double [][] Xp = new double[p+1][3];
		for (int idx=0; idx<X.length; idx++)
		{
			int n=0;
			Xp[0][0] = 1;
			Xp[0][1] = 1;
			Xp[0][2] = 1;
			for (int jdx=0; jdx<p; jdx++)
			{
				Xp[jdx+1][0] = Xp[jdx][0]*X[idx][0];
				Xp[jdx+1][1] = Xp[jdx][1]*X[idx][1];
				Xp[jdx+1][2] = Xp[jdx][2]*X[idx][2];
			} // for jdx
	
			for (int pdx=0; pdx <= p; pdx++)
			{
				for (int xp=pdx; xp >= 0; xp--)
				{
					for (int yp=pdx-xp; yp >= 0; yp--)
					{
						int zp = pdx-xp-yp;
						V[idx][n] = Xp[xp][0] * Xp[yp][1] * Xp[zp][2];
						n++;
					} // yp
				} // xp
			} // pdx
		} // idx
		return V;
	} // vander_3d

	public static final void derivative_2d(double[][] C_dx, double[][] C_dy, final double[][] C, final int n)
	{
		int s = 0;
		int m=0;
		for (int idx=1; idx<=n; idx++)
		{
			s = s+idx;
			for (int jdx=1; jdx<=idx; jdx++)
			{
				for (int k=0; k<C.length; k++)
				{
					C_dx[m][k] = (idx+1-jdx)*C[s+jdx-1][k];
					C_dy[m][k] = jdx*C[s+jdx][k];
				}
				m=m+1;
			} // for jdx
		} // for idx
	} // derivative_2d()

	// fist partial derivatives of the degree p polynomial basis functions (rows in C)
	public static final double[][] derivative_3d(final double [][] C, double [][] C_dx, double [][] C_dy, double [][] C_dz, final int p)
	{
		for (int idx=0; idx<C.length; idx++)
		{
			int nx=0;
			int ny=0;
			int nz=0;
			int n=0;
			for (int pdx=0; pdx<=p; pdx++)
			{
				for (int xp=pdx; xp>=0; xp--)
				{
					for (int yp=pdx-xp; yp>=0; yp--)
					{
						int zp=pdx-xp-yp;
						if (xp > 0)
						{
							C_dx[nx][idx] = xp*C[n][idx];
							nx++;
						}
						if (yp > 0)
						{
							C_dy[ny][idx] = yp*C[n][idx];
							ny++;
						}
						if (zp > 0)
						{
							C_dz[nz][idx] = zp*C[n][idx];
							nz++;
						}
						n++;
					} // yp
				} // xp
			} // pdx
		} // for idx
		return C_dx;
	} // derivative_3d

	public static final void display(final double A[][])
	{
		for (int i=0; i<A.length; i++)
		{
			for (int j=0; j<A[0].length; j++)
			{
        			System.out.print(A[i][j] + " ");
			} // for j
	        	System.out.println();
		} // for i
	}

	public static int [][] realloc2dInt(int [][] A, int l1, int l2)
	{
		int [][] B = new int[l1][l2];
		for (int i=0; i<Math.min(A[0].length,l2); i++)
		{
			//System.arraycopy(A[i], 0, B[i], 0, A.length);
			for (int j=0; j<Math.min(A.length,l1); j++)
			{
				B[j][i] = A[j][i];
			}
		}
		return B;
	}

	public static double [][] realloc2dDouble(double [][] A, int l1, int l2)
	{
		double [][] B = new double[l1][l2];
		for (int i=0; i<Math.min(A[0].length,l2); i++)
		{
			//System.arraycopy(A[i], 0, B[i], 0, A.length);
			for (int j=0; j<Math.min(A.length,l1); j++)
			{
				B[j][i] = A[j][i];
			}
		}
		return B;
	}

	// Wed Jul 25 16:54:58 MSK 2012
	// Karl Kästner, Berlin
	public static final double dist3(final double [] A, final double [] B)
	{
		double retval = Math.sqrt( 
			+ (B[0] - A[0])*(B[0] - A[0])
			+ (B[1] - A[1])*(B[1] - A[1])
			+ (B[2] - A[2])*(B[2] - A[2]) );
		return retval;
	} // dist3D

	public static final double dist2(final double [] A, final double [] B)
	{
		double retval = Math.sqrt( 
			+ (B[0] - A[0])*(B[0] - A[0])
			+ (B[1] - A[1])*(B[1] - A[1]) );
		return retval;
	} // dist3D

	// observe convergence of the sth-order derivative
	// see Strang, Fix Theorem 3.7, An Analysis of the Finite Element Method
	public static final int get_rate(final int m, final int k, final int s)
	{
		return Math.min(k-s, 2*(k-m));
		//double rate = 2*k;
	} // get_rate
} // class FEM

