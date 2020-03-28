// Mon Aug 13 21:17:48 MSK 2012
// Karl Kästner, Berlin

// benchmark for JAMA matrix inversion and multiplication
import Jama.*;

class Test_jama
{
	public static void time_inverse(int i_max, int j_max)
	{
		//int i_max 	= 5;
		//int j_max       = 100;

		int n = i_max;
		//for (int i=1; i<i_max; i++)
		{
			//int n = (2<<i);

			double [][] A  = new double[n][n];
			Matrix mA       = new Matrix(A);
			// create a test matrix
			identity(A);
			// invert the test matrix several times
			for (int j=0; j<j_max; j++)
			{
				mA.inverse();
			} // for j
		} // for i
	} // main

	public static int dummy(int i)
	{
		return i;
	}

	public static int dummy_loop(final int n)
	{
		int s = 0;
		for (int i=0; i<n; i++)
		{
			s += dummy(i);
		}
		return s;
	}

	public static void identity(double [][] A)
	{
		for (int i=0; i<A.length; i++)
		{
			for (int j=0; j<A[0].length; j++)
			{
				A[i][j] = 0.0;
			}
		}
		for (int i=0; i<Math.min(A.length,A[0].length); i++)
		{
			A[i][i] = 1;
		}
	}
} // class Test_jama

