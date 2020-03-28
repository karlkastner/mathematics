// Wed Jul 25 17:31:00 MSK 2012
// Karl Kästner, Berlin

public class Key2 extends Object
{
	public int [] A = new int[2];
	
	public Key2(int i0, int i1)
	{
		// sort values into the slots
		A[0] = Math.min(i0,i1);
		A[1] = i0+i1-A[0];
	}

	public boolean equals(Object OB)
	{
		if (!(OB instanceof Key2)) return false;
		Key2 KB = (Key2) OB;
		return ( A[0] == KB.A[0] && A[1] == KB.A[1]);
	} // hashCode

	public int hashCode()
	{
		return (A[0] + A[1]*65536);
	} // hashCode
} // Key3

