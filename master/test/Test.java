public class Test
{
	public static Object[] d3(double [][] C, double [][] dx, double [][] dy, double [][] dz, int n)
	{
		FEM.derivative_3d(C, dx, dy, dz, n);
		Object [] ret = new Object[3];
		ret[0] = dx;
		ret[1] = dy;
		ret[2] = dz;
		return ret;
	}
}

