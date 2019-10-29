

class Test_Object
{
	public static void main(String args[])
	{
//		int [] A = new int[3];
//		A[0] = 10; A[1] = 100; A[2] = 1000;
//		int [] B = new int[3];
//		B[0] = 10; B[1] = 100; B[2] = 1000;
//		int [] C = new int[3];
//		C[0] = 1; C[1] = 2; C[2] = 3;
		Key3 A = new Key3(10,100,1000);
		Key3 B = new Key3(1000,100,10);		
		Key3 C = new Key3(1,2,3);		
		
		Object oa = A;
		Object ob = B;
		Object oc = C;
		System.out.println(oa.equals(ob));
		System.out.println(oa.equals(oc));
		System.out.println(ob.equals(oc));
		System.out.println(oa.hashCode());
		System.out.println(ob.hashCode());
		System.out.println(oc.hashCode());
	}
}


