// Wed Jul 25 15:42:17 MSK 2012
// Karl Kästner, Berlin

public interface Tree
{
	public int get_vertex(int i, int j);
//	public int add_child(int p1, int p2); //number of arguments change in higher dimensions
	public int get_child(int i, int j);
	public boolean is_leave(int i);
	// splits a boundary element uniformly into several parts
	public int split(int i);
} // interface Tree

