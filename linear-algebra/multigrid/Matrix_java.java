import Jama.Matrix;
public class Matrix_java
{
	public static double det3(final double [][] A)
	{
		double det = (    A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
        			+ A[0][1]*(A[1][2]*A[2][0]  -A[1][0]*A[2][2])
		   		+ A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0])
				);
		return(det);
	}
	public static double [][] inv3(final double [][] A)
	{
	double det = det3(A);
	double [][] Ai = {{(A[1][1]*A[2][2] - A[1][2]*A[2][1])/det,(A[0][2]*A[2][1] - A[0][1]*A[2][2])/det,(A[0][1]*A[1][2] - A[0][2]*A[1][1])/det},
			  {	      (A[1][2]*A[2][0] - A[1][0]*A[2][2])/det,(A[0][0]*A[2][2] - A[0][2]*A[2][0])/det,(A[0][2]*A[1][0] - A[0][0]*A[1][2])/det},
			  { (A[1][0]*A[2][1] - A[1][1]*A[2][0])/det,(A[0][1]*A[2][0] - A[0][0]*A[2][1])/det,(A[0][0]*A[1][1] - A[0][1]*A[1][0])/det}};
	return(Ai);
	}

	public static double [] mldivide3(final double [][] A, final double [] b)
	{
	double det = det3(A);
	double [] x = {(A[1][1]*A[2][2] - A[1][2]*A[2][1])/det*b[0]+(A[0][2]*A[2][1] - A[0][1]*A[2][2])/det*b[1]+(A[0][1]*A[1][2] - A[0][2]*A[1][1])/det*b[2],
			  	      (A[1][2]*A[2][0] - A[1][0]*A[2][2])/det*b[0]+(A[0][0]*A[2][2] - A[0][2]*A[2][0])/det*b[1]+(A[0][2]*A[1][0] - A[0][0]*A[1][2])/det*b[2],
			   (A[1][0]*A[2][1] - A[1][1]*A[2][0])/det*b[0]+(A[0][1]*A[2][0] - A[0][0]*A[2][1])/det*b[1]+(A[0][0]*A[1][1] - A[0][1]*A[1][0])/det*b[2]};
	return(x);
	}

	// this is slow bc of cache inefficiency
	public static double [][] mtimes(final double [][] A, final double [][] B)
	{
		long startTime = System.currentTimeMillis();
		int n = A.length;
		double [][] C = new double [n][n];
		for (int i = 0; i<n;i++)
		{
		for (int j = 0; j<n; j++)
		{
		for (int k = 0; k<n; k++)
		{
			C[i][j] = C[i][j] + A[i][k]*B[k][j];
		}
		}
		}
		long endTime = System.currentTimeMillis();
		System.out.println("That took " + (endTime - startTime) + " milliseconds");
		return(C);
	}
	
	public static double [][] mtimes_jama(final double [][] A, final double [][] B)
	{
		Matrix mA = new Matrix(A);
		Matrix mB = new Matrix(B);
		Matrix mC = mA.times(mB);
		double [][] C = mC.getArray();
		return(C);
	}	

	public static double [][] inv_jama(final double [][] A)
	{
		Matrix matrix = new Matrix(A);
    		Matrix inverseMatrix = matrix.inverse();
	    	double[][] iA = inverseMatrix.getArray();
		return(iA);

	}
	public static double [][] mldivide_jama(final double [][] A, final double [][] b)
	{
		Matrix Am = new Matrix(A);
		Matrix bm = new Matrix(b);
		Matrix xm = Am.solve(bm);
	    	double[][] x = xm.getArray();
		return(x);
		
	}
	public static void do_nothing()
	{
	}
}

