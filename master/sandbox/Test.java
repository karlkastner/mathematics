// Sat Jun  2 23:00:20 MSK 2012
// Karl Kästner, Berlin
import java.io.*;
import java.util.*;

class Test
{
	public static void testdim()
	{
		double [][] A = new double[3][4];
		for (int i=0; i<A.length*A[0].length; i++)
		{
//			A[i] = (double) i;
		}
		/*
		for (int i=0; i<A.length; i++)
		{
			for (int j=0; j<A[0].length; j++)
			{
				A[i][j] = i*j;
			} // for j
		} // for i
		*/
		display(A);
	}

	public static double[][] increase(double A[][])
	{
		for (int i=0; i<A.length; i++)
		{
			for (int j=0; j<A[0].length; j++)
			{
        			A[i][j] += 1;
			} // for j
		} // for i
		return A;
	} // increase()

	public static void display(double A[][])
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

	public static void testincrease()
	{
		double A[][] = new double[2][3];
		for (int i=0; i<A.length; i++)
		{
			for (int j=0; j<A[0].length; j++)
			{
				A[i][j] = 0;
			} // for j
		} // for i

		display(A);
		increase(A);
		display(A);
	}



	public static void main(String[] args) throws InterruptedException
	{
		//testdim();
		//Test test = new Test();
		TestThread testThread = new TestThread();
		testThread.testThread(40000);
    	}
} // class Increase

