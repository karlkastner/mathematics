// Wed Jul 25 15:33:58 MSK 2012
// Karl Kästner, Berlin

import java.util.Hashtable;

// TODO, for 2D, also edges must be split
class Tree_2d  implements Tree
{
	private static final int VERTEX = 0; // 0:2
	private static final int EDGE   = 3; // 3:5
	private static final int NEIGHBOUR = 6; // 6:9
	private static final int CHILD  = 10; // 10:13
	private static final int LEN = 13;
	private int [][] data;
	private Hashtable<Key2,int[]> OpenHash;
	// number of elements
	private int n;
	Tree_1d E_tree;

	public Tree_2d()
	{
		n = 0;
		data = null;
	}

	public int get_vertex(int i, int j)
	{
		return data[i][VERTEX+j];
	} // get_vertex

	private void set_vertex(int i, int j, int val)
	{
		data[i][VERTEX+j] = val;
	} // set_vertex

	public int get_edge(int i, int j)
	{
		return data[i][EDGE+j];
	} // get_edge

	public void set_edge(int i, int j, int val)
	{
		data[i][EDGE+j] = val;
	} // set_edge

	public int get_neighbour(int i, int j)
	{
		return data[i][NEIGHBOUR+j];
	} // get_neighbour

	public void set_neighbour(int i, int j, int val)
	{
		data[i][NEIGHBOUR+j] = val;
	} // set_neighbour

	public int get_child(int i, int j)
	{
		return data[i][CHILD+j];
	} // get_child

	private void set_child(int i, int j, int val)
	{
		data[i][CHILD+j] = val;
	} // set_child
	
	public boolean is_leave(int i)
	{
		// is leave, if child is null
		return (0 == data[i][CHILD+0]);
	} // is_leave()
	
	double [][] get_all_points()
	{
		return E_tree.get_all_points();
	} // get_all_points

	public int get_all_leaves(int [][] T, int nt_[], int [][] B, int nb_[], double [][] P, int np_[], int [] R)
	{
		// reallocate memory
		// triangles
		T = new int[n][3];
		// boundary edges
		B = new int[n][2];
		// back reference
		R = new int[n];
		// number of triangles
		int nt = 0;
		// number of boundary edges
		int nb = 0;

		for (int i=0; i<n; i++)
		{
			// only generate if this is a leave
			if (is_leave(i))
			{
				// count splitted neighbours
				int ns = 0;
				for (int j=0; j<3; j++)
				{
					if (!E_tree.is_leave(get_edge(i,j))) ns++;
				} // for j
	
				int [] P_ = new int [3];
				P_[0] = get_vertex(i,0);
				P_[1] = get_vertex(i,1);
				P_[2] = get_vertex(i,2);
				switch (ns)
				{
					case 0 : // no edge is split
						// copy element vertices into the BC-array
						T[nt][0] = P_[0];
						T[nt][1] = P_[1];
						T[nt][2] = P_[2];
						nt++;
						// back association
						if (null != R) R[nt] = i;
						// edge opposit point 0	
						if (get_edge(i,0) < 0)
						{
							B[nb][0] = P_[1];
							B[nb][1] = P_[2];
							nb++;
						}
						// edge opposit point 1
						if (get_edge(i,0) < 0)
						{
							B[nb][0] = P_[0];
							B[nb][1] = P_[2];
							nb++;
						}
						// edge opposit point 2	
						if (get_edge(i,0) < 0)
						{
							B[nb][0] = P_[0];
							B[nb][1] = P_[1];
							nb++;
						}
					break;
					case 1 : // one edge was split
						// closure
						// connect hanging node with opposit side and write two triangles
	
						// rotate according to which side was split
						throw new RuntimeException("TODO");
						//break;
					case 2 : // two edges were split
						throw new RuntimeException("TODO");
						//break;
					default: // should not happen, as in this case the complete triangle
						// should have been split
						throw new RuntimeException();
				} // switch ns
			} // if is_leave
			// generate boundary elements
			throw new RuntimeException("TODO");
		} // for i
		// get point coordinates
		P = get_all_points();
		// write back values
		nt_[0] = nt;
		nb_[0] = nb;
		np_[0] = P.length;
		return 0;
	} // get_all_leaves()

	private int add_child(int p1, int p2, int p3)
	{
		// increase number of elements
		n=n+1;
		// set vertices
		set_vertex(n-1,0,p1);
		set_vertex(n-1,1,p1);
		set_vertex(n-1,2,p1);

		// Note: setting neighbours for 3D-boundaries is actually not necessary

		Key2 key;
		int[] val;
		// edge opposit point 0
		key = new Key2(p2, p3);
		val = OpenHash.get(key);
		if (null != val)
		{
			// set this boundary element as neighbour in 
			set_edge(val[0]-1, val[1], n);
			set_edge(n-1, 0, val[0]);
		} else {
			val = new int[2];
			val[0] = n;
			val[1] = 0;
			OpenHash.put(key, val);
		}

		// edge opposit point 1
		key = new Key2(p1, p3);
		val = OpenHash.get(key);
		if (null != val)
		{
			// set this boundary element as neighbour in 
			set_edge(val[0]-1, val[1], n);
			set_edge(n-1, 1, val[0]);
		} else {
			val = new int[2];
			val[0] = n;
			val[1] = 1;
			OpenHash.put(key, val);
		}

		// edge opposit point 2
		key = new Key2(p1, p2);
		val = OpenHash.get(key);
		if (null != val)
		{
			// set this boundary element as neighbour in 
			set_edge(val[0]-1, val[1], n);
			set_edge(n-1, 2, val[0]);
		} else {
			val = new int[2];
			val[0] = n;
			val[1] = 2;
			OpenHash.put(key, val);
		}

		// return 1-based index of new element
		return n;
	} // add_child

	// splits a triangular element uniformly into 4 parts
	public int split(int i)
	{
	    int retval = 0;		

	    // only split, this, if it is a leave
	    if (is_leave(i))
	    {
		int p1 = get_vertex(i,0);
		int p2 = get_vertex(i,1);
		int p3 = get_vertex(i,2);

		// split the three sides into halves and get the edge midpoints
		int p12 = E_tree.split(get_edge(i,0));
		int p13 = E_tree.split(get_edge(i,1));
		int p23 = E_tree.split(get_edge(i,2));

		// first child (exterior, first corner)
		set_child(i, 0 ,add_child( p1, p12, p13));
		// second child (exterior, second corner)
		set_child(i, 1, add_child(p12,  p2, p23));
		// third child (exterior, third corner)
		set_child(i, 2, add_child(p13, p23,  p3));
		// get the neighbour facing the boundary
		// fourth child (interior triangle)
		set_child(i, 2, add_child(p12, p23, p13));
		
		// recursively split the neighbours, if they have more than 2 split sides
		// Note: This is not required if this is a 3D boundary

		// test level of neighbour and recursively split the neighbours if necessary
		for (int j=0; j<3; j++)
		{
			// get neighbour
			int nj = get_neighbour(i,j);
			// test if this is a boundary
			if (nj < 0)
			{
				// Note: it is not required to explicitely store the boundaries
				// however, the children neighbouring the boundary have to get
				// their neigbourhood indices respectively
					// set_edge(val[0]-1, val[1], -child);
				throw new RuntimeException("TODO");

				/*
				
				
				// Note: not required for 3D boundary elements		

				// split the boundary triangle
				// (boundary sections are always split, if the neighbouring triangle were split)
				B_tree.split(-nj);
				// restore neighbourhood relations
				int child;
				Object val;
				Key2 key;
				// child 0
				child = B_tree.get_child(-nj,0);
				// get the neighbour facing the boundary
				key = new Key2(B_tree.get_vertex(child,0), B_tree.get_vertex(child,1));
				val = OpenHash.get(key);
				// set this boundary element as neighbour in 
				// boundary element has no back-pointer
				set_edge(val[0]-1, val[1], -child);
				// child 1
				child = B_tree.get_child(-nj,1);
				key = new Key2(B_tree.get_vertex(child,0), B_tree.get_vertex(child,1));
				val = OpenHash.get(key);
				set_edge(val[0]-1, val[1], -child);
				*/
			} else {
				// only split the neighbour if it was not yet split
				if (is_leave(nj))
				{
					// count number of split edges
					int ns = 0;
					for (int k=0; k<3; k++)
					{
						if (!E_tree.is_leave(get_edge(j,k))) ns++;
					} // for k
					
					// split neighbour, if it has more than 2 split edges
					if (ns > 2)
					{
						// recursive split
						split(nj);
					} // if ns > 3
				} // if not yet split
			} // if boundary
		} // for j (all neighbours)
		// successfully split
		retval = 1;
	    } // if is_leave
	    return retval;
	} // split
} // class Tree_2d

