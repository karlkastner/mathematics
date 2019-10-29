import java.lang.*;

class Sumbench
{
	public static void main(String [] argv)
	{
		// TODO, test with integer, repeat and estimate stdev instead of one long run (overflow)
		// TODO		axpy ( y = ax + y )
		float [][] A = new float[10000000][2];
		Rand(A); 
		long t;
		float sum;

		FArray2Ds FA = new FArray2Ds(A);

		t = java.lang.System.nanoTime();
		sum = sum_FA_man(FA);
		t = java.lang.System.nanoTime() - t;
		System.out.printf("%12.0f %12d\n",sum,t);

		t = java.lang.System.nanoTime();
		sum = sum_FA(FA);
		t = java.lang.System.nanoTime() - t;
		System.out.printf("%12.0f %12d\n",sum,t);

		t = java.lang.System.nanoTime();
		sum = sum2_man(A);
		t = java.lang.System.nanoTime() - t;
		System.out.printf("%12.0f %12d\n",sum,t);

		t = java.lang.System.nanoTime();
		sum = sum2(A);
		t = java.lang.System.nanoTime() - t;
		System.out.printf("%12.0f %12d\n",sum,t);

		t = java.lang.System.nanoTime();
		sum = sumn_fixed(A);
		t = java.lang.System.nanoTime() - t;
		System.out.printf("%12.0f %12d\n",sum,t);

		t = java.lang.System.nanoTime();
		sum = sumn_variable(A);
		t = java.lang.System.nanoTime() - t;
		System.out.printf("%12.0f %12d\n",sum,t);
	} // maind

	public static float sum_FA_man(FArray2Ds FA)
	{
		float sum = 0.0f;
		for (int i=0; i<FA.n; i++)
			sum += FA.at(i,0) + FA.at(i,1);
		return sum;
	}

	public static float sum_FA(FArray2Ds FA)
	{
		float sum = 0.0f;
		for (int i=0; i<FA.n; i++)
		for (int j=0; j<FA.m; j++)
			sum += FA.at(i,j);
		return sum;
	}

	public static float sum2_man(float [][] A) // [2][] generates error
	{
		float sum = 0.0f;
		for (int idx=0; idx<A.length; idx++)
		{
			sum += A[idx][0] + A[idx][1];
		}
		return sum;
	}

	public static float sum2(float [][] A) // [2][] generates error
	{
		float sum = 0.0f;
		for (int idx=0; idx<A.length; idx++)
		{
			for (int jdx=0; jdx<2; jdx++) // dangerous 
			{
				sum+=A[idx][jdx];
			}
		}
		return sum;
	}

	public static float sumn_fixed(float [][] A)
	{
		float sum = 0.0f;
		for (int idx=0; idx<A.length; idx++)
		for (int jdx=0; jdx<A[0].length; jdx++) // TODO dangerous
		{
			sum += A[idx][jdx];	
		}
		return sum;
	}

	public static float sumn_variable(float [][] A)
	{
		float sum = 0.0f;
		for (int idx=0; idx<A.length; idx++)
		for (int jdx=0; jdx<A[idx].length; jdx++)
		{
			sum += A[idx][jdx];	
		}
		return sum;
	}

	public static void Rand(float [][] A)
	{
		for (int i=0; i<A.length; i++)
		for (int j=0; j<A[i].length; j++)
			//A[i][j] = (float) Math.random();
			A[i][j] = 1; //(float) Math.random();
	}
} // Sumbench

// java to c translator

