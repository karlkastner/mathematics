public class stacktest {
 
	public static int recurse()
	{
		int ret = 0;
        try {
	    ret = recurse()+1;
	} catch (java.lang.StackOverflowError e) {
//	    System.out.print("Recursion depth on this system is " + i + ".");
	}
		return ret;
    }
 
    public static void main(String[] args) {
        recurse();
    }
}
 
