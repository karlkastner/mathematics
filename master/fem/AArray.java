// 2013-03-16 12:32:47.000000000 +0100
// Karl Kastner, Berlin

class AArray
{
/*
	public static int [][] resizeI(final int [][] T, final int n, final int m)
	{
		int [][] retval;
		int l;
		if (null == T)
		{
			n0 = 0;
		} else {
			n0 = T.length;
			// TODO, does this work if T[0] is empty?
			m0 = T[0].length;
		}

		if (l < n || )
		{
			int [][] T_new = new int[2*l][T[0].length];
			for (int i=0; i<l; i++)
			{
				for (int j=0; j<l; j++)
				{
					T_new[i][j] = T[i][j];
				} // j
			} // i
			retval =  T_new;
		} else {
			retval = T;
		}
		return retval;
	} // resizeI
*/

	// TODO this does not work for empy arrays, because second dimension is not known
	// so pass m as well
	public static int [][] resizeI(final int [][] T, final int n)
	{
		int [][] retval;
		int l;
		if (null == T)
		{
			l = 0;
		} else {
			l = T.length;
		}

		if (l < n)
		{
			int [][] T_new = new int[2*l][T[0].length];
			for (int i=0; i<l; i++)
			{
				for (int j=0; j<l; j++)
				{
					T_new[i][j] = T[i][j];
				} // j
			} // i
			retval =  T_new;
		} else {
			retval = T;
		}
		return retval;
	} // resizeI
	
	public static void swapI(int [] A, final int ia, final int ib)
	{
		int c = A[ia];
		A[ia] = A[ib];
		A[ib] = c;
	} // swapI
} // class AArray	

