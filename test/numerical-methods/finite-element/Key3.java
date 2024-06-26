// Wed Jul 25 13:28:21 MSK 2012
// Karl Kästner, Berlin

public class Key3 extends Object
{
	public int [] A = new int[3];
	
	public Key3(int i0, int i1, int i2)
	{
		// sort values into the slots
		A[2] = Math.min(i0,Math.min(i1,i2));
		A[0] = Math.max(i0,Math.max(i1,i2));
		A[1] = i0+i1+i2-A[0]-A[2];
	}

	public boolean equals(Object OB)
	{
		if (!(OB instanceof Key3)) return false;
		Key3 KB = (Key3) OB;
		return ( A[0] == KB.A[0] && A[1] == KB.A[1] && A[2] == KB.A[2]);
	} // hashCode

	public int hashCode()
	{
		return (A[0] + A[1]*256 + A[2]*65536);
	} // hashCode
} // Key3

