import java.util.concurrent.atomic.*;
import java.lang.ref.*;

class test_dist3d
{
	public static void main(String args[])
	{
		double [] A1 = {0.0,0,0};
		double [] B1 = {1,0,0};
		System.out.println(Dist.dist3d(A1,B1));
		double [] A2 = {0,1,0};
		double [] B2 = {0,0,0};
//		System.out.println(dist3d(A2,B2));
		double [] A3 = {0,0,1};
		double [] B3 = {0,0,0};
//		System.out.println(dist3d(A3,B3));

		Integer a = new Integer(0);
		inc(a);
		System.out.println(a);
	}

	public static void inc(Integer a)
	{
		a = a+1;
		System.out.println(a);
		// i = Integer.valueOf(i.intValue() + 1);  
	}
}

