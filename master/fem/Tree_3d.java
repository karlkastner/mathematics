// Wed Jul 25 11:56:35 MSK 2012
// Karl Kästner, Berlin

import java.util.Hashtable;
import java.util.Vector;
import java.util.Iterator;

// tree with non-conforming "red"-tetrahedra leaves and parents
class Tree_3d implements Tree
{
	public static final class INDEX
	{
		public final static int VERTEX = 0;
		public final static int EDGE = 4;
		public final static int NEIGHBOUR = 10;
		public final static int CHILD = 14;
		public final static int PARENT = 22;
		public final static int LEN = 23;
	}
/*	public static enum INDEX
	{
		VERTEX(0),
		EDGE(4),
		NEIGHBOUR(10),
		CHILD(14),
		PARENT(22),
		LEN(23);
		INDEX(int i)
		{
			index = i;
		};
		int index;
	}
*/

/*
	// data is stored in one array for efficiency reason
	// rows : triangle entries
	// col  0.. 3 : corner point indices
	private static final int VERTEX = 0; 
	// col  4.. 9 : edge indices
	private static final int EDGE   = 4;
        // col 10..13 : neighbours
	private static final int NEIGHBOUR = 10; 
	// col 14..21 : children
	private static final int CHILD  = 14;
	private static final int PARENT = 22;
	private final int LEN = 23;
*/
	private int [][] data;
	// number of elements including parents
	private int n;
	// number of elements before splitting
	private int n_;
	// back association from green triangles to non-conforming red parent
	private int [] R;
	// tree of edges
	// TODO : is a complete tree here required or can it be reduced like the boundary tree ?
	private Tree_1d E_tree;

	// hash of neighbours (faces)
	private Hashtable<Key3,int[]> N_hash;

	public Tree_3d(final double [][] P, final int [][] T, final int [][] Bc)
	{
		// allocate memory
		data   = new int[T.length][INDEX.LEN];
		N_hash = new Hashtable<Key3,int[]>();

		// add points
		E_tree = new Tree_1d(P);

		// back association array
		R = new int[T.length];

		// tetrahedra array to tree
		for (int i=0; i<T.length; i++)
		{
			// add tetrahedron
			add_child(T[i][0], T[i][1], T[i][2], T[i][3], 0, 0, 0, 0, 0);
			R[i] = i;
		} // for i (tetras)

		// resolve boundaries
		for (int i=0; i<Bc.length; i++)
		{
			Key3 key3 = new Key3(Bc[i][0], Bc[i][1], Bc[i][2]);
			int [] val = N_hash.get(key3);
			if (null != val)
			{
				// TOOD : instead of -1 there should be -1 the boundary group it belongs
				// to allow the application of different  boundary conditions
				set_neighbour(val[0]-1, val[1], -(i+1));
			} else {
				// should never happen
				//System.err.println(i + " " + Bc[i][0] + " " + Bc[i][1] + " " + Bc[i][2]);
				throw new RuntimeException();
			}
		} // for i (boundary triangles)

		/*
		if (!N_hash.isEmpty())
		{
			throw new RuntimeException("Tree_3d::Constructor: Inconsistent Mesh");
		} */
	} // constructor

	public final int get_vertex(final int i, final int j)
	{
		if ((DEBUG.level > 0) && ( j >= 4) ) throw new RuntimeException();
		// vertices have column indices 0-4
		return data[i][INDEX.VERTEX+j];
	} // get_vertex

	private final void set_vertex(final int i, final int j, final int val)
	{
		// vertices have column indices 0-4
		if ((DEBUG.level > 0) && ( j >= 4) ) throw new RuntimeException();
		data[i][INDEX.VERTEX+j] = val;
	} // set_vertex

	public final int get_edge(final int i, final int j)
	{
		if ((DEBUG.level > 0) && ( j >= 6) ) throw new RuntimeException();
		return data[i][INDEX.EDGE+j];
	} // get_edge

	private final void set_edge(final int i, final int j, final int val)
	{
		if ((DEBUG.level > 0) && ( j >= 6 )) throw new RuntimeException();
		data[i][INDEX.EDGE+j] = val;
	} // set_edge

	public final int get_neighbour(final int i, final int j)
	{
		if ((DEBUG.level > 0) && ( j >= 4) ) throw new RuntimeException();
		return data[i][INDEX.NEIGHBOUR+j];
	} // get_face

	public final void set_neighbour(final int i, final int j, final int val)
	{
		if ((DEBUG.level > 0) && ( j >= 4) ) throw new RuntimeException();
		data[i][INDEX.NEIGHBOUR+j] = val;
	} // set_face

	public final int get_child(final int i, final int j)
	{
		if ((DEBUG.level > 0) && ( j >= 8) ) throw new RuntimeException();
		return  data[i][INDEX.CHILD+j];
	} // get_child

	public final void set_child(final int i, final int j, final int val)
	{
		if ((DEBUG.level > 0) && ( j >= 8) ) throw new RuntimeException();
		data[i][INDEX.CHILD+j] = val;
	} // set_child

	public final int get_parent(final int i)
	{
		return data[i][INDEX.PARENT+0];
	} // get_parent

	public final void set_parent(final int i, final int val)
	{
		data[i][INDEX.PARENT+0] = val;
	} // set_parent

	public final boolean is_leave(final int i)
	{
		return (0 == data[i][INDEX.CHILD+0]);
	} // is_leave

	private final int add_child(final int p1, final int p2, final int p3, final int p4,
				final int n1, final int n2, final int n3, final int n4, final int parent)
	{
		// reallocate memory
		if (n+1 > data.length)
		{
			data = FEM.realloc2dInt(data, 2*data.length, INDEX.LEN);
		}

		// increase the number of elements
		n = n+1;
		// set the vertices
		set_vertex(n-1, 0, p1);
		set_vertex(n-1, 1, p2);
		set_vertex(n-1, 2, p3);
		set_vertex(n-1, 3, p4);
		set_parent(n-1, parent);
		Key2 key2;
		Key3 key3;
		int [] val;
		// resolve pointer to neighbour 0 (opposit point 0)
		key3 = new Key3(p2,p3,p4);
		val = N_hash.remove(key3);
		if (null != val)
		{
			set_neighbour(n-1, 0, val[0]);
			// only set back-link from the neighbour if it is not a boundary
			// TODO: actually neighbour will never be a boundary, superficial check
			if (val[0] > 0) set_neighbour(val[0]-1, val[1], n);
		} else {
			val = new int[3];
			val[0] =  n;
			val[1] =  0;
			val[2] = n1;
			N_hash.put(key3,val);
		}
		// resolve pointer to neighbour 1
		key3 = new Key3(p1,p3,p4);
		val  = N_hash.remove(key3);
		if (null != val)
		{
			set_neighbour(n-1, 1, val[0]);
			if (val[0] > 0) set_neighbour(val[0]-1, val[1], n);
		} else {
			val = new int[3];
			val[0] =  n;
			val[1] =  1;
			val[2] = n2;
			N_hash.put(key3,val);
		}
		// resolve pointer to neighbour 2
		key3 = new Key3(p1,p2,p4);
		val  = N_hash.remove(key3);
		if (null != val)
		{
			set_neighbour(n-1, 2, val[0]);
			if (val[0] > 0) set_neighbour(val[0]-1, val[1], n);
		} else {
			val = new int[3];
			val[0] =  n;
			val[1] =  2;
			val[2] = n3;
			N_hash.put(key3,val);
		}
		// resolve pointer to neighbour 3
		key3 = new Key3(p1,p2,p3);
		val  = N_hash.remove(key3);
		if (null != val)
		{
			set_neighbour(n-1, 3, val[0]);
			if (val[0] > 0)  set_neighbour(val[0]-1, val[1], n);
		} else {
			val = new int[3];
			val[0] =  n;
			val[1] =  3;
			val[2] = n4;
			N_hash.put(key3,val);
		}
		// resolve pointers to edges
		set_edge(n-1, 0, E_tree.resolve_edge(p1, p2, n) );
		set_edge(n-1, 1, E_tree.resolve_edge(p1, p3, n) );
		set_edge(n-1, 2, E_tree.resolve_edge(p1, p4, n) );
		set_edge(n-1, 3, E_tree.resolve_edge(p2, p3, n) );
		set_edge(n-1, 4, E_tree.resolve_edge(p2, p4, n) );
		set_edge(n-1, 5, E_tree.resolve_edge(p3, p4, n) );
		return n;
	} // add_child

	// gets all conforming leaves and green children of non-conforming leaves
	// generates the tetrahedron matrix P, T and Bc
	// to be used by the assembly routines
	// from the tetrahedron tree
	// includes closure routine
	public final Mesh_3d generate_mesh()
	{
		// TODO, count number of boundary elements
		// max nt is somewhat arbitrary, due to closure
		// allocate memory
		Mesh_3d mesh = new Mesh_3d( E_tree.np, 4*n, 4*n);
		// back association array
		R = new int[2*n];

		// copy point coordinates
		mesh.P = E_tree.get_all_points();
		mesh.np = E_tree.np;
		
		// for all nodes in the tetra tree
		for (int i=0; i<n; i++)
		{
			// if this is a leave
			if (is_leave(i))
			{
				// edge midpoint array
				int [] mp = new int[6];

				// count split edges 
				int s = 0;
				for (int j=0; j<6; j++)
				{
					int e = get_edge(i,j);
					if (!E_tree.is_leave(e-1))
					{
						s += (1<<j);
						mp[j] = E_tree.get_vertex(E_tree.get_child(e-1,0)-1,0);
					} else {
						mp[j] = 0;
					}
				} // for j

				// fetch point indices
				int p0 = get_vertex(i,0);
				int p1 = get_vertex(i,1);
				int p2 = get_vertex(i,2);
				int p3 = get_vertex(i,3);

				// fetch neighbour tetra indices
				int n0 = get_neighbour(i,0);
				int n1 = get_neighbour(i,1);
				int n2 = get_neighbour(i,2);
				int n3 = get_neighbour(i,3);

				// if no edge is split
				switch (s)
				{
					case 0: // conforming red triangle, write to output
						// add as red triangle to the T triangle-array
						mesh.add_T(p0,p1,p2,p3);
						// back association
						R[mesh.nt-1] = i;
						// generate boundary faces, if this element faces the boundary
						// face opposit point 0	
						if (n0 < 0) mesh.add_Bc(p1,p2,p3);
						// face opposit point 1
						if (n1 < 0) mesh.add_Bc(p0,p2,p3);
						// face opposit point 2	
						if (n2 < 0) mesh.add_Bc(p0,p1,p3);
						// face opposit point 3	
						if (n3 < 0) mesh.add_Bc(p0,p1,p2);
						//System.out.println("regular " + mesh.nt);
						break;
					case 1 :
						// closure : only edge 0 split (vertices 0-1)
						closure_2(mesh, p0, p1, p2, p3, mp[0], n0, n1, n2, n3, i);
						break;
					case 2 :
						// closure : only edge 1 split (vertices 0-2)
						closure_2(mesh, p0, p2, p1, p3, mp[1], n0, n2, n1, n3, i);
						break;
					case 4 : 
						// closure : only edge 2 split (vertices 0-3)
						closure_2(mesh, p0, p3, p2, p1, mp[2], n0, n3, n2, n1, i);
						break;
					case 8 :
						// closure : only edge 3 split (vertices 1-2)
						closure_2(mesh, p1, p2, p0, p3, mp[3], n1, n2, n0, n3, i);
						break;
					case 16 :
						// closure : only edge 4 split (vertices 1-3)
						closure_2(mesh, p1, p3, p2, p0, mp[4], n1, n3, n2, n0, i);
						break;
					case 32 :
						// closure : only edge 5 split (vertices 2-3)
						closure_2(mesh, p2, p3, p0, p1, mp[5], n2, n3, n0, n1, i);
						break;
					case 11 : // 2^0 + 2^1 + 2^3
						// closure, face 3 split (opposit point 3) (0 based)
						closure_4(mesh, p0, p1, p2, p3, mp[0], mp[1], mp[3], n0, n1, n2, n3, i);
						break;
					case 21 : // 2^0 + 2^2 + 2^4
						// closure, face 2 split (opposit point 2) (0 based)
						closure_4(mesh, p0, p1, p3, p2, mp[0], mp[2], mp[4], n0, n1, n3, n2, i);
						break;
					case 38 : // 2^1 + 2^2 + 2^5
						// closure, face 1 split (opposit point 1) (0 based)
						closure_4(mesh, p0, p3, p2, p1, mp[2], mp[1], mp[5], n0, n3, n2, n1, i);
						break;
					case 56 : // p12-p13-p23, e3-e4-e5, 2^5 + 2^4 + 2^3
						// closure, face 0 split (opposit point 0) (0 based)
						closure_4(mesh, p3, p1, p2, p0, mp[4], mp[5], mp[3], n3, n1, n2, n0, i);
						break;
					case 3 : // p01-p02, e0-e1 1 + 2 = 3
						closure_3(mesh, p0, p1, p2, p3, mp[0], mp[1], n0, n1, n2, n3, i);
						break;
					case 5 : // p01-p03, e0-e2 1 + 4 = 5
						closure_3(mesh, p0, p1, p3, p2, mp[0], mp[2], n0, n1, n3, n2, i);
						break;
					case  9 : // p01-p12, e0-e3 1 + 8 = 9
						closure_3(mesh, p1, p0, p2, p3, mp[0], mp[3], n1, n0, n2, n3, i);
						break;
					case 17 : // p01-p13, e0-e4 1 + 16 =17
						closure_3(mesh, p1, p0, p3, p2, mp[0], mp[4], n1, n0, n3, n2, i);
						break;
					//                    e0-e5, not neighbouring
					case  6 : // p02-p03, e1-e2, 2 + 4 = 6
						closure_3(mesh, p0, p2, p3, p1, mp[1], mp[2], n0, n2, n3, n1, i);
						break;
					case 10 : // p02-p12, e1-e3,  2 + 8 = 10
						closure_3(mesh, p2, p0, p1, p3, mp[1], mp[3], n2, n0, n1, n3, i);
						break;
					//      p02-p13, e1-e4, 2 + 16 = 18	(not neighbouring)
					case 34 : // p02-p23, e1-e5, 2 + 32 = 34	
						closure_3(mesh, p2, p0, p3, p1, mp[1], mp[5], n2, n0, n3, n1, i);
						break;
					//               e2-e3, not neighbouring
					case 20 : // p03-p13, e2-e4, 4 + 16 =20
						closure_3(mesh, p3, p0, p1, p2, mp[2], mp[4], n3, n0, n1, n2, i);
						break;
					case 36 : // p03-p23, e2-e5, 4 + 32 =36
						closure_3(mesh, p3, p0, p2, p1, mp[2], mp[5], n3, n0, n2, n1, i);
						break;
					case 24 : // p12-p13, e3-e4, 8 + 16  =24
						closure_3(mesh, p1, p2, p3, p0, mp[3], mp[4], n1, n2, n3, n0, i);
						break;
					case 40 : // p12-p23, e3-e5, 8 + 32 =40
						closure_3(mesh, p2, p1, p3, p0, mp[3], mp[5], n2, n1, n3, n0, i);
						break;
					case 48 : // p13-p23, e4-e5, 16 + 32 =48
						closure_3(mesh, p3, p1, p2, p0, mp[4], mp[5], n3, n1, n2, n0, i);
						break;
					default:
						// three shoud not happen, as then the triangle should be split regularly
						// more than three should not happpen, as there are only three edges
						throw new RuntimeException(i + " " + s + " " +Integer.toBinaryString(s));
				} // if no edge split yet
			} // if not a leave
		} // for int i (tetrahedra)
		return mesh;
	} // get_all_leves


	// splits a tetrahedron into 8 children
	// recursively splits the neighbours, if necessary
	public final int split(final int i)
	{
	    int retval = 0;

	    // check if child is zero
	    if (is_leave(i))
	    {
		// mark this element as split
		set_child(i, 0, Integer.MAX_VALUE);
		
		// conformity rule 1: make sure, that the neighbours of the parent are split as well
		int parent_ = get_parent(i);
		if (parent_ > 0)
		{
			for (int j=0; j<4; j++)
			{
				int nj = get_neighbour(parent_-1, j);
				if (0 == nj) throw new RuntimeException();
				if (nj > 0)
					split(nj-1);
			}
		} // parent > 0
/*
		for (int j=0; j<4; j++)
		{ 
			// check if there is not yet a neighbour assigned to that side

			if (DEBUG.level > 0)
			{
				final int [][] q = {
					{ 1, 2, 3},
					{ 0, 2, 3},
					{ 0, 1, 3},
					{ 0, 1, 2} };
				int p0 = get_vertex(i,q[j][0]);
				int p1 = get_vertex(i,q[j][1]);
				int p2 = get_vertex(i,q[j][2]);
				Key3 key3 = new Key3(p0, p1, p2);
				int [] val = N_hash.get(key3);
				if ( ( null == val ) && ( 0 == get_neighbour(i,j) ) )
				{
					// does never happen
					throw new RuntimeException();
				}
				if ( ( 0 == get_neighbour(i,j) ) )
				{
					if (0 == get_parent(i))
					{
						// does never happen
						throw new RuntimeException();
					}
					int nj = get_neighbour(get_parent(i)-1,j);
					if (nj != val[2])
					{
						throw new RuntimeException(nj + " " + val[2]);
					}
				}
			} // DEBUG checks

			if (0 == get_neighbour(i,j))
			{
				// there is not yet a neighbour, so recursively split its parent
				int parent = get_parent(i);
				if (parent > 0)
				{
					// to be neighbour is a tetrahedron
					int nj = get_neighbour(parent-1,j);
					split(nj - 1);
				} else {
					// to be neighbour is a boundary - that should never happen
					// as neighbouring boundaries are split emmideately
					throw new RuntimeException();
				}
			} // null != val
		} // for j (conformity rule 1)
*/

		// conformity rule 2:
		// at first split the unsplit elements adjacent to the edge to be split
		for (int j=0; j<6; j++)
		{
			int parent = E_tree.get_parent(get_edge(i,j)-1);
			if (0 != parent)
			{
				//System.out.println(parent-1);
				Vector<Integer> v = E_tree.get_adjacent_element_v(parent-1);
				int help [] = new int[v.size()];
				Iterator<Integer> it = v.iterator();
				int k =0;
				while (it.hasNext())
				{
					int nj = (int)it.next();
					help[k] = nj;
					k++;
				}
				// avoid concurrent array modifiaction
				for (k=0; k<help.length; k++)
				{
					if (help[k] > 0) split(help[k]-1);
				}
			}
		}

		// get the point indices
		int p1 = get_vertex(i,0);
	 	int p2 = get_vertex(i,1);
	        int p3 = get_vertex(i,2);
	        int p4 = get_vertex(i,3);

		// get neighbours
		int n1 = get_neighbour(i,0);
		int n2 = get_neighbour(i,1);
		int n3 = get_neighbour(i,2);
		int n4 = get_neighbour(i,3);

		// split the six edges if not yet split, else get the midpoints
		int p12 = E_tree.split(get_edge(i,0)-1); // p1, p2
		int p13 = E_tree.split(get_edge(i,1)-1); // p1, p3
		int p14 = E_tree.split(get_edge(i,2)-1); // p1, p4
		int p23 = E_tree.split(get_edge(i,3)-1); // p2, p3
		int p24 = E_tree.split(get_edge(i,4)-1); // p2, p4
		int p34 = E_tree.split(get_edge(i,5)-1); // p3, p4

		// children 0:3, exterior tetrahedra
		// neighbours of inner edges are immediately resolved
		// no need to put correct neighbours into hash values
		// can be reduced to hand ofver parent id, if ordering is specified as here
//		set_child(i, 0, add_child( p1, p12, p13, p14, 0, n2, n3, n4));
//		set_child(i, 1, add_child( p2, p12, p23, p24, 0, n1, n3, n4));
//		set_child(i, 2, add_child( p3, p13, p23, p34, 0, n1, n3, n4));
//		set_child(i, 3, add_child( p4, p14, p24, p34, 0, n1, n2, n3));
		set_child(i, 0, add_child(  p1, p12, p13, p14,  0, n2, n3, n4, i+1));
		set_child(i, 1, add_child( p12,  p2, p23, p24, n1,  0, n3, n4, i+1));
		set_child(i, 2, add_child( p13, p23,  p3, p34, n1, n2,  0, n4, i+1));
		set_child(i, 3, add_child( p14, p24, p34,  p4, n1, n2, n3,  0, i+1));

	        // the four interior tetrahedra (three possibilities)
		// choose shortest edge as common edge
		double d1234 = FEM.dist3(E_tree.get_point(p12-1), E_tree.get_point(p34-1));
		double d1324 = FEM.dist3(E_tree.get_point(p13-1), E_tree.get_point(p24-1));
		double d1423 = FEM.dist3(E_tree.get_point(p14-1), E_tree.get_point(p23-1));
		// TODO, pass neighbours correctly
		if (d1423 < d1324 && d1423 < d1234)
		{
//	            set_child(i, 4, add_child(p12, p13, p14, p23, 0, 0, 0, 0)); // 4
	            set_child(i, 4, add_child(p14, p23, p24, p34, n1,  0,  0,  0, i+1)); // 1
		    set_child(i, 5, add_child(p13, p23, p14, p34,  0, n2,  0,  0, i+1)); // 2
		    set_child(i, 6, add_child(p12, p14, p23, p24,  0,  0, n3,  0, i+1)); // 3
	            set_child(i, 7, add_child(p12, p13, p23, p14,  0,  0,  0, n4, i+1)); // 4
		} else if (d1234 < d1324) {
//		    set_child(i, 4, add_child(p12, p34, p13, p23,  0, 0, 0, 0));
//		    set_child(i, 5, add_child(p12, p34, p13, p14,  0, 0, 0, 0)); // 2
//		    set_child(i, 6, add_child(p12, p34, p14, p24,  0, 0, 0, 0)); // 
//		    set_child(i, 7, add_child(p12, p34, p23, p24, n1, 0, 0, 0)); // 1
		    set_child(i, 4, add_child(p12, p23, p24, p34, n1,  0,  0,  0, i+1)); // 1
		    set_child(i, 5, add_child(p13, p12, p14, p34,  0, n2,  0,  0, i+1)); // 2
		    set_child(i, 6, add_child(p12, p14, p34, p24,  0,  0, n3,  0, i+1)); // 3
		    set_child(i, 7, add_child(p12, p13, p23, p34,  0,  0,  0, n4, i+1)); // 4
		} else {
//		    set_child(i, 4, add_child(p13, p24, p23, p34, 0, 0, 0, 0));
//		    set_child(i, 5, add_child(p13, p24, p14, p34, 0, 0, 0, 0));
//		    set_child(i, 6, add_child(p13, p24, p12, p23, 0, 0, 0, 0));
//		    set_child(i, 7, add_child(p13, p24, p12, p14, 0, 0, 0, 0));
		    set_child(i, 4, add_child(p13, p23, p24, p34, n1,  0,  0,  0, i+1)); // 1
		    set_child(i, 5, add_child(p13, p24, p14, p34,  0, n2,  0,  0, i+1)); // 2
		    set_child(i, 6, add_child(p12, p14, p13, p24,  0,  0, n3,  0, i+1)); // 3
		    set_child(i, 7, add_child(p12, p13, p23, p24,  0,  0,  0, n4, i+1)); // 4
		} // shortest edge

		// test level of neighbour and recursively split the neighbours if necessary
		for (int j=0; j<4; j++)
		{
			// get neighbour
			int nj = get_neighbour(i,j);
			// check if neighbour is a boundary and split
			if (nj < 0)
			{
				// neighbour is a boundary
				switch (j)
				{
					case 0: 
						split_boundary(p2, p3, p4, p23, p24, p34, nj);
					break;
					case 1:
						split_boundary(p1, p3, p4, p13, p14, p34, nj);
					break;
					case 2:
						split_boundary(p1, p2, p4, p12, p14, p24, nj);
					break;
					case 3:
						split_boundary(p1, p2, p3, p12, p13, p23, nj);
					break;
					default:
						// should not happen
						throw new RuntimeException();
				} // switch j
			} // if nj < 0
		} // for j (neighbours via faces)
		
		// neighbours via edges
		for (int j=0; j<6; j++)
		{
			int e = get_edge(i,j);
			Vector<Integer> v = E_tree.get_adjacent_element_v(e-1);
			int help [] = new int[v.size()];
			int l = 0;
			Iterator<Integer> it = v.iterator();
			while (it.hasNext())
			{
				help[l] = 0;
				int nj = (int)it.next();
				// check if the neighbour has to be splitted
				// only split the neighbour if it was not yet split
				if (is_leave(nj-1))
				{
					// count number of split edges
					int s = 0;
					for (int k=0; k<6; k++)
					{
						if (!E_tree.is_leave(get_edge(nj-1,k)-1)) s += (1<<k);
					} // for k
					// tetrahedra with more than three split edges are split as well
					// neigbours differ at most by 1 level
					// sort out cases with closure routine
					switch (s)
					{
						case  0 : // no closure
						// single edge
						case  1 : // edge 0
						case  2 : // edge 1
						case  4 : // edge 2
						case  8 : // edge 3
						case 16 : // edge 4
						case 32 : // edge 5
						// three edges on one face
						case 11 : // face 3 (edge 0, 1, 3) 1+2+ 8
						case 21 : // face 2 (edge 0, 2, 4) 1+4+16
						case 38 : // face 1 (edge 1, 2, 5) 2+4+32
						case 56 : // face 0 (edge 3, 4, 5) 8+16+32
						// two neighbouring edges
						case  3 : // p01-p02, e0-e1,  1 +  2 =  3
						case  5 : // p01-p03, e0-e2,  1 +  4 =  5
						case  9 : // p01-p12, e0-e3,  1 +  8 =  9
						case 17 : // p01-p13, e0-e4,  1 + 16 = 17
						case  6 : // p02-p03, e1-e2,  2 +  4 =  6
						case 10 : // p02-p11, e1-e3,  2 +  8 = 10
						case 34 : // p02-p23, e1-e5,  2 + 32 = 34
						case 20 : // p03-p13, e2-e4,  4 + 16 = 20
						case 36 : // p03-p23, e2-e5,  4 + 32 = 36
						case 24 : // p12-p13, e3-e4,  8 + 16 = 24
						case 40 : // p12-p23, e3-e5,  8 + 32 = 40
						case 48 : // p13-p23, e4-e5, 16 + 32 = 48
						//System.out.println("tetra " + i + " n " + (nj-1) + " " + s  + " Not splitting (side)");
						break;
						default:
						////System.out.println("Splitting (side): " + i + " " + (nj-1) + " " + s);
						//System.out.println("tetra " + i + " n " + (nj-1) + " " + s  + " Splitting (side)");
						help[l] = nj;
					} // if ns > 3

				} // if not yet split
				l++;
			} // while has next
			// avoid concurrent array modifiaction run time error
			for (int k=0; k<help.length; k++)
			{
				if (help[k] > 0) split(help[k]-1);
			}
		} // for j

		// successfully split, return 1
                retval = 1;
	    } // is a leave
 	    // TODO, make this a void function
	    return retval;
	} // split

	private final void split_boundary(final int p1, final int p2, final int p3,
				final int p12, final int p13, final int p23, final int nj)
	{
		// resolve neighbour pointers of the children facing the the domain boundary
		Key3 key3;
		int [] val;
		// resolve pointer to neighbour 0
		key3 = new Key3(p1,p12,p13);
		val = N_hash.remove(key3);
		if (null != val)
		{
			// mark neihgbour in child as a boundary
			set_neighbour(val[0]-1, val[1], nj);
		} else {
			// should not happen
			throw new RuntimeException();
		}
		// resolve pointer to neighbour 1
		key3 = new Key3(p2,p12,p23);
		val  = N_hash.remove(key3);
		if (null != val)
		{
			set_neighbour(val[0]-1, val[1], nj);
		} else {
			throw new RuntimeException();
		}
		// resolve pointer to neighbour 2
		key3 = new Key3(p3,p13,p23);
		val  = N_hash.remove(key3);
		if (null != val)
		{
			set_neighbour(val[0]-1, val[1], nj);
		} else {
			throw new RuntimeException();
		}
		// resolve pointer to neighbour 3
		key3 = new Key3(p12,p13,p23);
		val  = N_hash.remove(key3);
		if (null != val)
		{
			set_neighbour(val[0]-1, val[1], nj);
		} else {
			throw new RuntimeException();
		}
	} // split_boundary()

	// splits marked triangles and neighbours
	// M : marked elements for refinement, 1 based indices
	public final void refine(final int [] M)
	{
		n_ = n;
		// split each marked cell
		try {
		for (int i=0; i<M.length; i++)
		{
			// translate red/green T-entry into red-parent Tree-entry
			int j = R[M[i]-1];
			// split the triangles in the Tree and recursively its neighbours if required
			split(j);
		} // for i
		} catch (java.lang.StackOverflowError e) {
			e.printStackTrace();
			System.err.println(e);
			throw new RuntimeException(" " + n);
		}
		if ((DEBUG.level > 0)) verify();
		if ((DEBUG.level > 0)) check_neighbour();
	} // refine()

	//
	public final void coarsen(final int [] M)
	{
		// coarsening: select triangle with lowest error, remove the three points, remesh the region
		// do not remove a triangle where the neighbour was removed (colour the mesh rgb)
	} // coarsen

	public final void verify()
	{
		boolean change;
		do
		{
			change = false;
		for (int i=0; i<n; i++)
		{
			if (is_leave(i))
			{
			// count number of split edges
			int s = 0;
			for (int k=0; k<6; k++)
			{
				if (!E_tree.is_leave(get_edge(i,k)-1)) s += (1<<k);
			} // for k
			// tetrahedra with more than three split edges are split as well
			// neigbours differ at most by 1 level
			// sort out cases with closure routine
			switch (s)
			{
				case  0 : // no closure
				// a single edge
				case  1 : // edge 0
				case  2 : // edge 1
				case  4 : // edge 2
				case  8 : // edge 3
				case 16 : // edge 4
				case 32 : // edge 5
				// three edges enclosing one face split
				case 11 : // face 3 (edge 0, 1, 3) 1+2+ 8
				case 21 : // face 2 (edge 0, 2, 4) 1+4+16
				case 38 : // face 1 (edge 1, 2, 5) 2+4+32
				case 56 : // face 0 (edge 3, 4, 5) 8+16+32
				// two neighbouring edges
				case  3 : // p01-p02, e0-e1,  1 +  2 =  3
				case  5 : // p01-p03, e0-e2,  1 +  4 =  5
				case  9 : // p01-p12, e0-e3,  1 +  8 =  9
				case 17 : // p01-p13, e0-e4,  1 + 16 = 17
				//           p01-p23, e0-35   1 + 32 = 33 (not neighbouring)
				case  6 : // p02-p03, e2-e2,  2 +  4 =  6
				case 10 : // p02-p12, e1-e3,  2 +  8 = 10
				//           p02-p13, e1-e4,  2 + 16 = 18 (not neighbouring)
				case 34 : // p02-p23, e1-e5,  2 + 32 = 34
				//           p03-p12, e2-e3,  4 +  8 = 12 (not neighbouring)
				case 20 : // p03-p13, e2-e4,  4 + 16 = 20
				case 36 : // p03-p23, e2-e5,  4 + 32 = 36
				case 24 : // p12-p13, e3-e4,  8 + 16 = 24
				case 40 : // p12-p23, e3-e5,  8 + 32 = 40
				case 48 : // p13-p23, e4-e5, 16 + 32 = 48
				//System.out.println("tetra " + i + " n " + (nj-1) + " " + s  + " Not splitting (side)");
				break;
				default:
				////System.out.println("Splitting (side): " + i + " " + (nj-1) + " " + s);
				//System.out.println("tetra " + i + " n " + (nj-1) + " " + s  + " Splitting (side)");
				//System.err.println("inconsistency: "
				//	+ i + " " + s + " " + Integer.toBinaryString(s));
				//split(i);
				change = true;
				throw new RuntimeException("inconsistency: "
					+ i + " " + s + " " + Integer.toBinaryString(s));
				//split(i);
			} // if ns > 3
			} // if not yet split
		} // for i
		} while (change);
	} // verify

	// covers the six cases where one edge is split
	// exemplified for edge 0, p0-p1
	public final void closure_2(Mesh_3d mesh,
		final int p0, final int p1, final int p2, final int p3,
		final int p01,
		final int n0, final int n1, final int n2, final int n3, final int i)
	{
		if ((DEBUG.level > 0) && (0 == p0*p1*p2*p3*p01)) throw new RuntimeException();

		// tetrahedra
		mesh.add_T( p0, p01, p2, p3); R[mesh.nt-1] = i;
		mesh.add_T(p01,  p1, p2, p3); R[mesh.nt-1] = i;
//		System.out.println("split 2 " + mesh.nt + " " + (mesh.nt-1));

		// boundaries
		// face opposit point 0	
		if (n0 < 0)
		{
			mesh.add_Bc(p1,p2,p3);
		}
		// face opposit point 1
		if (n1 < 0)
		{
			mesh.add_Bc(p0,p2,p3);
		}
		// face opposit point 2	
		if (n2 < 0)
		{
			mesh.add_Bc(p0,p01,p3);
			mesh.add_Bc(p1,p01,p3);
		}
		// face opposit point 3	
		if (n3 < 0)
		{
			mesh.add_Bc(p1,p2,p01);
			mesh.add_Bc(p0,p2,p01);
		}
	} // closure_2 

	// green irregular refinement
	// splits an tetrahedra into three tetrahedra
	// this covers the 12 cases where two neighbouring edges are split
	// not intended for the 3 cases where two edges in opposition are split
	// there are two possible tesselations
	private final void closure_3(Mesh_3d mesh,
		final int p0, final int p1, final int p2, final int p3,
		final int p01, final int p02,
		final int n0, final int n1, final int n2, final int n3, final int i)
	{
		if ((DEBUG.level > 0) && (0 == p0*p1*p2*p3*p01*p02)) throw new RuntimeException();

		if (p01 < p02)
		{
			// child tetrahedra
			mesh.add_T(p0, p01, p02, p3); R[mesh.nt-1] = i;
			mesh.add_T(p1, p01, p02, p3); R[mesh.nt-1] = i;
			mesh.add_T(p1,  p2, p02, p3); R[mesh.nt-1] = i;
//			System.out.println("split 3 " + mesh.nt + " " + (mesh.nt-1) + " " + (mesh.nt-2));
			// face zero opposit point zero is boundary
			if (n0 < 0)
			{
				mesh.add_Bc(p1, p2, p3);
			}
			// face one opposit point one is a boundary
			if (n1 < 0)
			{
				mesh.add_Bc(p0, p02, p3);
				mesh.add_Bc(p2, p02, p3);
			}
			// face two opposit point 2 is a boundary
			if (n2 < 0)
			{
				mesh.add_Bc(p0, p01, p3);
				mesh.add_Bc(p1, p01, p3);
			}
			// face three opposit point 3 is a boundary
			if (n3 < 0)
			{
				mesh.add_Bc(p0, p01, p02);
				mesh.add_Bc(p1, p01, p02);
				mesh.add_Bc(p1,  p2, p02);
			}
		} else {
			// swap p01 with p02 and p1 with p2 an n1 with n2
			closure_3(mesh, p0, p2, p1, p3, p02, p01, n0, n2, n1, n3, i);
		}
	} // closure_3

	// covers the four cases where three edges on one face are split
	// not the four cases where three edges are split which are not edges of the same face
	// exemple cases covers splitting of face 3, opposit point 3
	public final void closure_4(Mesh_3d mesh,
		final int p0, final int p1, final int p2, final int p3,
		final int p01, final int p02, final int p12,
		final int n0, final int n1, final int n2, final int n3, final int i)
	{
		if ((DEBUG.level > 0) && (0 == p0) || (0 == p1) || (0 == p2) || (0 == p3) || (0 == p01) || (0 == p02) || (0 == p12))
			throw new RuntimeException(p0 + " " + p1 + " " + p2 + " " + p3 + " " + p01 + " " + p02 + " " + p12);

		// tetrahedra
		mesh.add_T( p0, p01, p02, p3); R[mesh.nt-1] = i;
		mesh.add_T(p01,  p1, p12, p3); R[mesh.nt-1] = i;
		mesh.add_T(p02, p12,  p2, p3); R[mesh.nt-1] = i;
		mesh.add_T(p01, p12, p02, p3); R[mesh.nt-1] = i;
//		System.out.println("split 4 " + mesh.nt + " " + (mesh.nt-1) + " " + (mesh.nt-2) + " " + (mesh.nt-3));
		
		// boundaries
		if (n0 < 0)
		{
			mesh.add_Bc(p3,p1,p12);
			mesh.add_Bc(p3,p2,p12);
		}
		if (n1 < 0)
		{
			mesh.add_Bc(p3,p0,p02);
			mesh.add_Bc(p3,p2,p02);
		}
		if (n2 < 0)
		{
			mesh.add_Bc(p3,p0,p01);
			mesh.add_Bc(p3,p1,p01);
		}
		// yes it is unlikely, but indeed possible that the three split face is a boundary
		if (n3 < 0)
		{
			mesh.add_Bc(p01, p02, p12);
			mesh.add_Bc( p0, p01, p12);
			mesh.add_Bc( p1, p01, p12);
			mesh.add_Bc( p2, p02, p12);
		}
	} // closure_4

	// test function
	public final void check_neighbour()
	{
		boolean err = false;
		Hashtable<Key3,Integer> H = new Hashtable<Key3,Integer>();
		// verifies neighbourhood links
		for (int i=0; i<n; i++)
		{
			Key3 key;	
			Integer val;
			// side 0
			key = new Key3(	get_vertex(i,1), get_vertex(i,2), get_vertex(i,3));
			val = H.get(key);
			if ( null == val )
			{
				H.put(key,i+1);
			} else {
				if (val != get_neighbour(i,0))
				{
					System.err.println(i+1 + " " + val + " " + get_neighbour(i,0));
					err = true;
				}
			}
			// side 1
			key = new Key3(	get_vertex(i,0), get_vertex(i,2), get_vertex(i,3));
			val = H.get(key);
			if ( null == val )
			{
				H.put(key,i+1);
			} else {
				if (val != get_neighbour(i,1))
				{
					System.err.println(i+1 + " " + val + " " + get_neighbour(i,1));
					err = true;
				}
			}
			// side 2
			key = new Key3(	get_vertex(i,0), get_vertex(i,1), get_vertex(i,3));
			val = H.get(key);
			if ( null == val )
			{
				H.put(key,i+1);
			} else {
				if (val != get_neighbour(i,2))
				{
					System.err.println(i + " " + val + " " + get_neighbour(i,2));
					err = true;
				}
			}
			// side 3
			key = new Key3(	get_vertex(i,0), get_vertex(i,1), get_vertex(i,2));
			val = H.get(key);
			if ( null == val )
			{
				H.put(key,i+1);
			} else {
				if (val != get_neighbour(i,3))
				{
					System.err.println(i + " " + val + " " + get_neighbour(i,3));
					err = true;
				}
			}
		}
		if (err) throw new RuntimeException();
	} // check_neighbour()
} // class Tree_3D

