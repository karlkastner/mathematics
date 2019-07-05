// Mon 12 Dec 21:00:49 CET 2016
// resample a function
class Resample
{
	// double [] last;
	//double [] id;
	//public static int [] resample_d_min(final double [] x, final double dmin)
	//public static int [] resample_d_min(final double [] x, final double d_min)
	public static int [] resample_d_min(final double [] x, final double d_min)
	{
		// allocate memory
		int [] id = new int[x.length];
		id[0]=1;
		int last=0;
		for (int i=1; i<x.length; i++)
		{
			id[i] = id[i-1];
			if (Math.abs(x[i]-x[last]) >= d_min)
			{
				last  = i;
				id[i] = id[i]+1;
			}
		} // for i	
		return id;
	} // resample_d_min
} // class Resample

