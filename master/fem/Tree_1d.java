// Wed Jul 25 15:16:13 MSK 2012
// Karl Kästner, Berlin

import java.util.Hashtable;
import java.util.Vector;

public final class Tree_1d implements Tree
{
	// column 0 : first vertex
	// column 1 : second vertex
	// column 2 : first child
	// column 3 : second child
	private static final int VERTEX = 0;
	private static final int CHILD  = 2;
	private static final int PARENT = 4;
	private static final int LEN = 5;
	private static final int DIM = 3; // TODO, what about 2D ?
	private int [][] data;
	private int n;
	private Hashtable<Key2,int[]> hash;

	// point coordinates
	private double P[][];
	// number of points
	public int np; // TODO privatise
	// elements being adjacent to edges
	private Object [] adjacent_element_v;

	public Tree_1d( final double [][] P )
	{
		np = P.length;
		// allocate memory
		hash = new Hashtable<Key2,int[]>();
		this.P = new double[np][P[0].length]; // TODO dangerous
		data = new int[1][LEN];
		
		n = 0;
		for (int i=0; i<np; i++)
		{
			// TODO, use ArrayCopy
			for (int j=0; j<P[i].length; j++)
			{
				this.P[i][j] = P[i][j];
			}
		}
		adjacent_element_v = new Object[np];
	} // Tree_1d

	public final int get_vertex(final int i, final int j)
	{
		return data[i][VERTEX+j];
	} // get_vertex

	private final void set_vertex(final int i, final int j, final int val)
	{
		data[i][VERTEX+j] = val;
	} // set_vertex

	private final void set_child(final int i, final int j, final int val)
	{
		data[i][CHILD+j] = val;
	} // set_child

	public final int get_child(final int i, final int j)
	{
		return data[i][CHILD+j];
	} // get_child

	public final Vector<Integer> get_adjacent_element_v(final int i)
	{
		return (Vector<Integer>) adjacent_element_v[i];
	} // get_adjacent_element_v
	
	public final void add_adjacent_element(final int i, final int element)
	{
		((Vector<Integer>)adjacent_element_v[i]).add(new Integer(element));	
	} // add_adjacent_element

	public final int get_parent(final int i)
	{
		return data[i][PARENT+0];
	} // get_parent

	private final void set_parent(final int i, final int val)
	{
		data[i][PARENT+0] = val;
	} // set_parent

	// p1, p2 and paretnt are expected 1-based
	public final int add_child(final int p1, final int p2, final int parent)
	{
		// check if it exists (only necessary in combination with 3d)
		Key2 key2 = new Key2(p1,p2);
		int[] val = hash.get(key2);
		int retval;
		
		if (null == val)
		{
			// reallocate memory
			if (n+1 > data.length)
			{
				data = FEM.realloc2dInt(data, 2*data.length, LEN);
				Object[] ha = new Object[2*adjacent_element_v.length];
				for (int i=0; i<n; i++)
				{
					ha[i] = adjacent_element_v[i];
				}
				adjacent_element_v = ha;
			}
	
			set_vertex(n, 0, p1);
			set_vertex(n, 1, p2);
			set_parent(n, parent);
			adjacent_element_v[n] = new Vector<Integer>();
			n++;
			// add to hash
			val = new int[1];
			val[0] = n;
			hash.put(key2, val);
			retval = n;
		} else {
			// already exists
			retval = val[0];
		}
		return retval;
	} // add_child

	// adds segment midpoint
	// expects 1-based indices
	private final int add_point(final int p1, final int p2)
	{
		// reallocate memory if necessary
		if (np+1 > P.length)
		{
			P = FEM.realloc2dDouble(P, 2*P.length, DIM);
		} // if reallocation necessary

		// same routine for 1, 2 and 3D
		for (int i=0; i<P[np].length; i++)
		{
			P[np][i] = 0.5*(P[p1-1][i]+P[p2-1][i]);
		}
		np++;
		return np;
	} // add_point

	public final double [] get_point(final int i)
	{
		return P[i];
	} // get_point

	// get point coordinates
	public final double [][] get_all_points()
	{
		return P;
	} // get_all_points

	public final boolean is_leave(final int i)
	{
		// element is a leave, if the children are null
		return (0 == data[i][CHILD+0]);
	} // is_leave

	// p1 and p2 are expected 1-based
	public final int resolve_edge(final int p1, final int p2, final int element)
	{
		// edges cannot be removed from the hash,
		// as they may belong to more than two elements
		Key2 key2 = new Key2(p1, p2);
		int [] val  = hash.get(key2);
		int edge;
		if (null != val)
		{
			// edge does already exist
 			edge = val[0];
		} else {
			// usually edges where allready created during splitting
			// however, this is not so at load, in this case, the parent is zero
			// and for interior edges of the faces 
//			System.out.println("adding edge");
			edge = add_child(p1, p2, 0);
		}
		add_adjacent_element(edge-1, element);
		return edge;
	} // resolv edge

	// splits the edge, if not yet split and returns the midpoint index
	public final int split(final int i)
	{
		int pm;
		if (!is_leave(i))
		{
			// edge already split, midpoint exists already
			// midpoint is per definition vertex 0 in the children
			pm = get_vertex(get_child(i,0)-1,0);
		} else {
			int p1 = get_vertex(i, 0);
			int p2 = get_vertex(i, 1);
			
			// create a new mid-point
			pm = add_point(p1, p2);

			// split edge segment (add_child expects 1-based indices)
			int e1 = add_child(pm, p1, i+1);
			// children are stored with 1-based indices
			set_child(i, 0, e1);
			int e2 = add_child(pm, p2, i+1);
			set_child(i, 1, e2);
		} // if !is_leave
		return pm;
	} // split_2
} // class Tree_1d

