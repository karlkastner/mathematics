// 2014-03-15 19:28:48.193847503 +0700
// Karl Kastner, Berlin

class Swap
{
	public static void swap( int[] a, final int i, final int j) //int [] b)
	{
		//int tmp = a[0];
		//b[0] = a[0];
		//a[0] = tmp;
		int tmp = a[i];
		a[i] = a[j];
		a[j] = tmp;
	} // swap
} // Swap

