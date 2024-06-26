//Fri Jul 27 18:54:02 MSK 2012
// Karl Kästner, Berlin

public final class Potential_3D_Coulomb implements Potential_3D
{
		public final double[] func(final double[][] q)
		{
			double[] f = new double[q.length];
			for (int i=0; i<q.length; i++)
			{
				f[i] = -1.0/Math.sqrt( q[i][0]*q[i][0] + q[i][1]*q[i][1] + q[i][2]*q[i][2] );
			}
			return f;
		} // f
} // class Func_Coulomb


