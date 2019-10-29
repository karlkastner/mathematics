public class FArray2D<T> // actually fortran style
{
	private
	T [] array;
	int n;
	int m;

	java.util.Hashtable ha = new java.util.Hashtable();

	public FArray2D(final int n, final int m) //, Class<T> type)
	{
		this.n = n;
		this.m = m;
//		Class c = (Class) ((java.lang.reflect.ParameterizedType).getClass().getGenericSuperclass());
		Class c = (Class<T>) ((java.lang.reflect.ParameterizedType) getClass()
                            .getGenericSuperclass()).getActualTypeArguments()[0];
		array = (T[]) java.lang.reflect.Array.newInstance(c, m*n); //arrayType.getComponentType( ), list.size( )));  
		//array = new T[n*m];
	} // FArray2D

	// TODO dangerous, assumes, that each row in A is of equal length
	public FArray2D(T [][] A)
	{
		m = A[0].length;
		n = A.length;
		//array = new T[n*m];
		Class c = (Class<T>) ((java.lang.reflect.ParameterizedType) getClass()
                            .getGenericSuperclass()).getActualTypeArguments()[0];
		array = (T[]) java.lang.reflect.Array.newInstance(c, m*n);
		for (int y=0; y<n; y++)
		{
			for (int x=0; x<m; x++)
			{
				array[y*m+x] = A[y][x];
			} // x
		} // y
	} // FArray2D
	
	public T at(final int x, final int y)
	{
		return array[x+y*m];
	}
} // class FArray2D

