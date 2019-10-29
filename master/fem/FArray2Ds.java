public class FArray2Ds // actually fortran style
{
	private
	float [] array;
	final int n;
	final int m;

	java.util.Hashtable ha = new java.util.Hashtable();

	public FArray2Ds(final int n, final int m) //, Class<T> type)
	{
		this.n = n;
		this.m = m;
		array = new float[n*m];
	} // FArray2D

	// TODO dangerous, assumes, that each row in A is of equal length
	public FArray2Ds(float [][] A)
	{
		m = A[0].length;
		n = A.length;
		array = new float[n*m];
		for (int y=0; y<n; y++)
		{
			for (int x=0; x<m; x++)
			{
				array[y*m+x] = A[y][x];
			} // x
		} // y
	} // FArray2D
	
	public float at(final int y, final int x)
	{
		return array[x+y*m];
	}
} // class FArray2D

