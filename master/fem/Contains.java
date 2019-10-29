// 2013-03-16 23:01:26 +0700
// Karl Kastner, Berlin

import Jama.*;

public class Contains
{



	// returns coefficients xc, yc, r^2
	public static double [] circle2D(final double [][] P, final int [][] T, final int tdx)
	{
		double x0 = P[T[tdx][0]][0];
		double y0 = P[T[tdx][0]][1];
		double x1 = P[T[tdx][1]][0];
		double y1 = P[T[tdx][1]][1];
		double x2 = P[T[tdx][2]][0];
		double y2 = P[T[tdx][2]][1];
		final double [][] A = { { 2*x0, 2*y0, 1 },
		  		    { 2*x1, 2*y1, 1 },
                                    { 2*x2, 2*y2, 1 }};
		final double [][] b = { { x0*x0 + y0*y0 }, { x1*x1 + y1*y1 }, { x2*x2 + y2*y2 } };
		
		Matrix mA = new Matrix(A);
		Matrix mB = new Matrix(b);
		
		Matrix mC = mA.solve(mB);

		double [][] C = mC.getArray();
		double [] retval = new double[3];
		retval[0] = C[0][0];
		retval[1] = C[1][0];
		retval[2] = -C[2][0] + retval[0]*retval[0] + retval[1]*retval[1];

		return retval;
	} // circle

	// TODO, dummy function for compilation
	public static double [] circle2D(double [][] P, double [][] T, int tdx)
	{
		return new double [0];
	}

	// circle of circumference
	public static boolean conains2d_circ(double [][] P, int [][] T, int tdx, int pdx)
	{
		double [] circle = circle2D(P, T, tdx);
		double dx = P[pdx][0] - circle[0];
		double dy = P[pdx][1] - circle[1];

		boolean retval = (dx*dx + dy*dy) < circle[2];

		return retval;
	}

	public static void main (String [] args)//test()
	{
		int [][] T_2d = {{0, 1, 2}};
		double [][] P_2d = {{0,0}, {0,1}, {1,0}, {0.5,0.5}, {0.25,0.25}, {0.75,0.75}};

/*		
		System.out.println(
					contains2D(P_2d,T_2d,0,3) + " " +
					contains2D(P_2d,T_2d,0,4) + " " +
					contains2D(P_2d,T_2d,0,5) );
*/

		int [][] T_3d = {{0, 1, 2, 3}};
		double [][] P_3d = {{0,0,0}, {0,0,1}, {0,1,0}, {1,0,0}, {0.5,0.5,0.0}, {0.25,0.25,0.25}, {0.75,0.75,0.75}};

/*		
		System.out.println(
					contains3D(P_3d,T_3d,0,4) + " " +
					contains3D(P_3d,T_3d,0,5) + " " +
					contains3D(P_3d,T_3d,0,6) );
*/
                                   
		double [] circle = circle2D(P_2d, T_2d, 0);

		System.out.println(circle[0] + " " + circle[1] + " " + circle[2]);
	} // test
} // Contains

