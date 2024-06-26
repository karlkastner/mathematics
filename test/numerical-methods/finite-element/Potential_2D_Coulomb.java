// Sun Jun  3 16:27:56 MSK 2012
// Karl Kästner, Berlin
//import Func.*;

public class Potential_2D_Coulomb implements Potential_2D
{
		public double[] func(double[][] q)
		{
			double[] f = new double[q.length];
			for (int i=0; i<q.length; i++)
			{
				f[i] = -1.0/Math.sqrt( q[i][0]*q[i][0] + q[i][1]*q[i][1]);
			}
			return f;
		} // f
} // class Func_Coulomb


