	public static boolean contains2D(final double [][] P, final int [][] T,
							final int tdx, final int pdx)
	{
		double A[][] = { { 1.0, 1.0, 1.0 },
       	                  { P[T[tdx][0]][0], P[T[tdx][1]][0], P[T[tdx][2]][0] },
       	                  { P[T[tdx][0]][1], P[T[tdx][1]][1], P[T[tdx][2]][1] } };
		double [][] pv = {{1.0}, {P[pdx][0]}, {P[pdx][1]}};
		Matrix mA    = new Matrix(A);
		Matrix  mP    = new Matrix(pv);

		Matrix mC = mA.solve(mP);
		double [][] c = mC.getArray();	

		boolean retval = (c[0][0] >= 0.0) && (c[1][0] >= 0.0) && (c[1][0] >= 0.0);

		return retval;	
	} // contains2D

