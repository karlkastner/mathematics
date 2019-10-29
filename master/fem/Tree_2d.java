// Wed Jul 25 15:33:58 MSK 2012
// Karl Kästner, Berlin

import java.util.Hashtable;
import java.util.Arrays;

// TODO triangle tree can also be build from arbitrary triangulations, but then triangles overlap
//	(minimum triangle that contains all points)
//	simple, but not optimal: get encircling circle and make this circle incircle of a triangle


// TODO, for 2D, also edges must be split
public final class Tree_2d  implements Tree
{
	public static final class INDEX
	{
		public static final int VERTEX    = 0;  // 0:2   (3)
		public static final int EDGE      = 3;  // 3:5   (3)
		public static final int NEIGHBOUR = 6;  // 6:8   (3)
		public static final int CHILD     = 9;  // 9:12 (4) 
		public static final int PARENT    = 13; // 13    (1)
		public static final int LEN       = 14; // total
	}
	private int [][] data;
	// number of elements
	private int n;
	// number of elements before splitting
	private int n_;
	// tree of edges
	Tree_1d E_tree;

	// back association from green triangles to non-conforming red parent
	private int [] R;

	// hash of neighbours (edges)
	//private Hashtable<Key2,int[]> OpenHash;
	private Hashtable<Key2,int[]> N_hash;

	// maximum edges that can be split before the element is split itself
	// 0 : causes uniform refinement
	// 1 : only 1-split auxiliary triangles
	// 2 : 1 and 2-split auxiliary triangles
	public int MAX_SPLIT_EDGES = 1;

	public Tree_2d()
	{
	}

	public Tree_2d(final double [][] P, final int [][] T, final int [][] Bc)
	{
		// allocate memory
		data   = new int[T.length][INDEX.LEN];
		N_hash = new Hashtable<Key2,int[]>();

		// add points
		E_tree = new Tree_1d(P);

		// back association array
		R = new int[T.length];

		// triangle array to tree
		for (int i=0; i<T.length; i++)
		{
			// add triangle
			add_child(T[i][0], T[i][1], T[i][2], 0, 0, 0, 0); //, 0, 0, 0, 0, 0);
 			R[i] = i;
		} // for i (tetras)

		// resolve boundaries
		for (int i=0; i<Bc.length; i++)
		{
			Key2  key2 = new Key2(Bc[i][0], Bc[i][1]);
			int [] val = N_hash.get(key2);
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

	public final int [][] get_data()
	{
		return data;
	}

	public final Tree_1d get_E_tree()
	{
		return E_tree;
	}	

	public final int get_vertex(final int i, final int j)
	{
		return data[i][INDEX.VERTEX+j];
	} // get_vertex

	private void set_vertex(final int i, final int j, final int val)
	{
		data[i][INDEX.VERTEX+j] = val;
	} // set_vertex

	public final int get_edge(final int i, final int j)
	{
		return data[i][INDEX.EDGE+j];
	} // get_edge

	public void set_edge(final int i, final int j, final int val)
	{
		data[i][INDEX.EDGE+j] = val;
	} // set_edge

	public final int get_neighbour(final int i, final int j)
	{
		return data[i][INDEX.NEIGHBOUR+j];
	} // get_neighbour

	public void set_neighbour(final int i, final int j, final int val)
	{
		data[i][INDEX.NEIGHBOUR+j] = val;
	} // set_neighbour

	public final int get_child(final int i, final int j)
	{
		// TODO improve by data[i][CHILD] + j;
		return data[i][INDEX.CHILD+j];
	} // get_child

	private void set_child(final int i, final int j, final int val)
	{
		data[i][INDEX.CHILD+j] = val;
	} // set_child

	public final int get_parent(final int i)
	{
		return data[i][INDEX.PARENT+0];
	} // get_parent

	public final void set_parent(final int i, final int j)
	{
		data[i][INDEX.PARENT+0] = j;
	} // get_parent
	
	public final boolean is_leaf(final int i)
	{
		// is leaf, if child is null
		return (0 == data[i][INDEX.CHILD+0]);
	} // is_leaf()
	
	final double [][] get_all_points()
	{
		return E_tree.get_all_points();
	} // get_all_points

	public final int [] get_R()
	{
		return R;
	} // get_all_points

	// was get_all_leafs
	public final Mesh_2d generate_mesh(boolean close)
	{
		// int [][] T, int nt_[], int [][] B, int nb_[], double [][] P, int np_[], int [] R)
		Mesh_2d mesh = new Mesh_2d(E_tree.np, 4*n, 4*n);

		// back association array
		R = new int[2*n];

		// copy point coordinates
		mesh.P  = E_tree.get_all_points();
		mesh.np = E_tree.np;
/*
		// reallocate memory
		// triangles
		int [][] T = new int[n][3];
		// boundary edges
		B = new int[n][2];
		// back reference
		R = new int[n];
		// number of triangles
		int nt = 0;
		// number of boundary edges
		int nb = 0;
*/
		for (int i=0; i<n; i++)
		{
			// only generate if this is a leaf
			if (is_leaf(i))
			{
				// fetch point indices
				int p0 = get_vertex(i,0);
				int p1 = get_vertex(i,1);
				int p2 = get_vertex(i,2);

				if (!close)
				{
					mesh.add_T(p0,p1,p2);
				} else {

				// count splitted neighbours
				int mp[] = new int[3];
				int ns = 0;
				for (int j=0; j<3; j++)
				{
					int e = get_edge(i,j);
					if (!E_tree.is_leaf(e-1)) 
					{
						// TODO edges only have one child, last zero is spurious
						int ec = E_tree.get_child(e-1,0);
//						System.out.println((i+1) + " " + j + " " + e + " " + ec);
						mp[j] = E_tree.get_vertex(ec-1,0);
						ns  += 1<<j;
//						System.out.println(j);
					}
				} // for j
//				System.out.println("ns " + ns + " " + Arrays.toString(mp));



				// fetch neighbour triangle indices
				int n0 = get_neighbour(i,0);
				int n1 = get_neighbour(i,1);
				int n2 = get_neighbour(i,2);

				// side opposit point 0	
				if (n0 < 0) mesh.add_Bc(p1,p2);
				// side opposit point 1
				if (n1 < 0) mesh.add_Bc(p0,p2);
				// side opposit point 2	
				if (n2 < 0) mesh.add_Bc(p0,p1);
	
				//int [] P_ = new int [3];
				//P_[0] = get_vertex(i,0);
				//P_[1] = get_vertex(i,1);
				//P_[2] = get_vertex(i,2);
				switch (ns)
				{
					case 0 : // no edge is split
					         // conforming red triangle, write to output
						 // add as red triangle to the T triangle-array
						 mesh.add_T(p0,p1,p2);

						// back association
						R[mesh.nt-1] = i;

						// generate boundary edges, if this element faces the boundary


/*

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
*/
					break;
					// 1-edge split closures
					// connect hanging node with opposit side and write two triangles
					case 1 : // only edge 1 is split
						mesh.add_T(p0,p1,mp[0]);
						R[mesh.nt-1] = i;
						mesh.add_T(p0,mp[0],p2);
						R[mesh.nt-1] = i;
						break;
					case 2 : // only edge 2 is split
						mesh.add_T(p0,p1,mp[1]);
						R[mesh.nt-1] = i;
						mesh.add_T(p1,p2,mp[1]);
						R[mesh.nt-1] = i;
						break;
					case 4 : // only edge 3 was split
						mesh.add_T(p0,mp[2],p2);	
						R[mesh.nt-1] = i;
						mesh.add_T(p1,p2,mp[2]);	
						R[mesh.nt-1] = i;
						break;
					// two edge split closure
					case 3 :
						throw new RuntimeException("TODO 3");
						// break;	
					case 5 :
						throw new RuntimeException("TODO 5");
						// break;	
					case 6 :
//						System.out.println("Closing " + (i+1));
						throw new RuntimeException("TODO 6");
						//break;
					default: // should not happen, as in this case the complete triangle
						// should have been split
						throw new RuntimeException("TODO" + ns);
				} // switch ns
				} // if close
			} // if is_leaf
			// TODO automatically generate neighbours
		} // for i
		// get point coordinates
		//P = get_all_points();
		// write back values
		//nt_[0] = nt;
		//nb_[0] = nb;
		//np_[0] = P.length;
		return mesh;
	} // generate_mesh()

	// TODO, n of parrent does not need to be passed as an argument
	private final int add_child(final int p1, final int p2, final int p3, 
                    final int n1, final int n2, final int n3, final int parent)
	{
		// reallocate memory
		if (n+1 > data.length)
		{
			data = FEM.realloc2dInt(data, 2*data.length, INDEX.LEN);
		}

		// increase number of elements
		n=n+1;
	
		// set vertices
		set_vertex(n-1,0,p1);
		set_vertex(n-1,1,p2);
		set_vertex(n-1,2,p3);
	
		// set parent
		set_parent(n-1, parent);

		Key2 key;
		int[] val;

		// edge opposit p1
		key = new Key2(p2, p3);
		val = N_hash.get(key);
		if (null != val)
		{
			set_neighbour(     n-1,      0, val[0]);
			set_neighbour(val[0]-1, val[1],      n);

//			set_edge(val[0]-1, val[1],      n);
//			set_edge(     n-1,      0, val[0]);
		} else {
			// add open edge to hash
			val = new int[2];
			val[0] = n;
			val[1] = 0;
			N_hash.put(key, val);
		}

		// edge opposit p2
		key = new Key2(p1, p3);
		val = N_hash.get(key);
		if (null != val)
		{
			set_neighbour(     n-1,      1, val[0]);
			set_neighbour(val[0]-1, val[1],      n);
				

//			set_edge(val[0]-1, val[1],      n);
//			set_edge(     n-1,      1, val[0]);
		} else {
			val = new int[2];
			val[0] = n;
			val[1] = 1;
			N_hash.put(key, val);
		}

		// edge opposit p3
		key = new Key2(p1, p2);
		val = N_hash.get(key);
		if (null != val)
		{
			set_neighbour(     n-1,      2, val[0]);
			set_neighbour(val[0]-1, val[1],      n);

//			set_edge(val[0]-1, val[1], n);
//			set_edge(n-1, 2, val[0]);
		} else {
			val = new int[2];
			val[0] = n;
			val[1] = 2;
			N_hash.put(key, val);
		}
		// resolve pointers to edges
		// TODO this is not explicitely necessary in 2d, if this
		// handled separately during initial mesh set up
		set_edge(n-1, 0, E_tree.resolve_edge(p2, p3, n) );
		set_edge(n-1, 1, E_tree.resolve_edge(p1, p3, n) );
		set_edge(n-1, 2, E_tree.resolve_edge(p1, p2, n) );

		// return 1-based index of new element
		return n;
	} // add_child

	// splits marked triangles and neighbours
	// M : marked elements for refinement, 1 based indices
	public final void refine(final int [] M)
	{
		n_ = n;
		// split each marked element
		try {
		for (int i=0; i<M.length; i++)
		{
			// translate red/green T-entry into red-parent Tree-entry
			int j = R[M[i]-1];
//			System.out.println("Requesting Splitting triangle " + (j+1) + "(" + M[i] + ")");
			// split the triangles in the Tree and recursively its neighbours if required
			split(j);
		} // for i
		} catch (java.lang.StackOverflowError e) {
			e.printStackTrace();
			System.err.println(e);
			throw new RuntimeException(" " + n);
		}
//		if ((DEBUG.level > 0)) verify();
//		if ((DEBUG.level > 0)) check_neighbour();
	} // refine()

	// splits a triangular element uniformly into 4 parts
	public int split(int i)
	{
	    int retval = 0;		

	    // only split, this, if it is a leaf
	    if (is_leaf(i))
	    {
//		System.out.println("Splitting triangle " + (i+1));
//%		System.out.println(Arrays.deepAsList(get_data()));
//		System.out.println(Arrays.deepToString(get_data()));
//		ArrayPrint2d(get_data());
		// mark this element temporarily as split to avoid infinite recursion
		set_child(i, 0, Integer.MAX_VALUE);

		// conformity rule 1 : make sure, that the neighbours of the parent are split as well
		int parent_ = get_parent(i);
		if (parent_ > 0)
		{
			for (int j=0; j<3; j++)
			{
				int nj = get_neighbour(parent_-1, j);
//				System.out.println(i + " " + parent_ + " " + j + " " + nj);
				if (0 == nj) throw new RuntimeException((i+1) + " " + (parent_)  + " " + j + " Zero neighbour");
				if (nj > 0)
					split(nj-1);
			}
		} // parent > 0

		// conformity rule 2 does not exist in 2d

		int p1 = get_vertex(i,0);
		int p2 = get_vertex(i,1);
		int p3 = get_vertex(i,2);

		int n1 = get_neighbour(i,0);
		int n2 = get_neighbour(i,1);
		int n3 = get_neighbour(i,2);

		// split the three sides in halve and get the edge midpoints
		int e1  = get_edge(i,0);
		int p23 = E_tree.split(e1-1);
		int e2  = get_edge(i,1);
		int p13 = E_tree.split(e2-1);
		int e3  = get_edge(i,2);
//		System.out.println(e3 + " " + i + " " + E_tree.get_n() );
		int p12 = E_tree.split(e3-1);

		// first child (exterior, first corner)
		int c1 = add_child( p1, p12, p13, 0, n2, n3, i+1);
		set_child(i, 0 ,c1);
		// second child (exterior, second corner)
		//int c2 = add_child(p12,  p2, p23, n1, 0, n3, i+1);
		int c2 = add_child(p2, p23, p12, n1, 0, n3, i+1);
		set_child(i, 1, c2);
		// third child (exterior, third corner)
		//int c3 =  add_child(p13, p23,  p3, n1, n2, 0, i+1);
		int c3 =  add_child(p3, p13, p23, n1, n2, 0, i+1);
		set_child(i, 2, c3);
		// get the neighbour facing the boundary
		// fourth child (interior triangle)
		int c4 = add_child(p12, p23, p13, 0, 0, 0, i+1);
		set_child(i, 3, c4);
	
		// split boundaries
		for (int j=0; j<3; j++)
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
						split_boundary(p2, p3, p23, nj);
						//split_boundary(p1, p2, p12, nj);
					break;
					case 1:
						split_boundary(p1, p3, p13, nj);
						//split_boundary(p1, p3, p13, nj);
					break;
					case 2:
						split_boundary(p1, p2, p12, nj);
						//split_boundary(p2, p3, p23, nj);
					break;
					default:
						// should not happen
						throw new RuntimeException();
				} // switch j
			} // if nj < 0
		} // for j (neighbours via edges)
	
		// recursively split the neighbours, if they have more than 2 split sides
		// Note: This is not required if this is a 3D boundary

		// test level of neighbour and recursively split the neighbours if necessary
		for (int j=0; j<3; j++)
		{
			// get neighbour
			int nj = get_neighbour(i,j);
//			System.out.println("T " + (i+1) + " neighbour " + nj);
			// test if this is a boundary
			if (nj < 0)
			{
				// if the neighbour is a boundary, nothing has to be done

/*
				// Note: it is not required to explicitely store the boundaries
				// however, the children neighbouring the boundary have to get
				// their neigbourhood indices respectively
					// set_edge(val[0]-1, val[1], -child);
				// throw new RuntimeException("TODO");
				
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
				val = N_hash.get(key);
				// set this boundary element as neighbour in 
				// boundary element has no back-pointer
				set_edge(val[0]-1, val[1], -child);
				// child 1
				child = B_tree.get_child(-nj,1);
				key = new Key2(B_tree.get_vertex(child,0), B_tree.get_vertex(child,1));
				val = N_hash.get(key);
				set_edge(val[0]-1, val[1], -child);
*/
			} else if (nj > 0 ){
				// only split the neighbour if it was not yet split
				if (is_leaf(nj-1))
				{
					// count number of split edges
					int ns = 0;
					for (int k=0; k<3; k++)
					{
						// TODO in 2d there is a short cut, as edges link uniquely the neighbours of the triangle
						// neighbours can be checked with is leaf 
						int e_ = get_edge(nj-1,k);
						if (!E_tree.is_leaf(e_-1)) ns++;
						//if (!E_tree.is_leaf(get_edge(j,k))) ns++;
					} // for k
					
					// split neighbour, if it has more than MAX_SPLIT edges
					// this has either to be one or two
//					System.out.println(ns);
					if (ns > MAX_SPLIT_EDGES)
					{
						// recursive split
						split(nj-1);
					} // if ns > 3
				} // if not yet split
			} // if not a boundary
		} // for j (all neighbours)

		// successfully split
		retval = 1;
	    } // if is_leaf
	    return retval;
	} // split

	private final void split_boundary(final int p1, final int p2,
				final int p12, final int nj)
	{
		// resolve neighbour pointers of the children facing the the domain boundary
		Key2 key;
		int [] val;
		// resolve pointer to neighbour 0
		key = new Key2(p1,p12);
		val = N_hash.remove(key);
		if (null != val)
		{
			// mark neihgbour in child as a boundary
			set_neighbour(val[0]-1, val[1], nj);
		} else {
			// should not happen
			throw new RuntimeException();
		}
		// resolve pointer to neighbour 1
		key  = new Key2(p2,p12);
		val  = N_hash.remove(key);
		if (null != val)
		{
			set_neighbour(val[0]-1, val[1], nj);
		} else {
			throw new RuntimeException();
		}
	} // split_boundary()

	public void refine_nonobtuse(final double abstol)
	{
		int n0 = n;
		// first loop, refine all obtuse triangles where the longest edge is on the boundary
		// TODO can be sped up by only testing the boundary triangles
		for (int i=0; i<n0; i++)
		{
			if (is_obtuse(i,abstol))
			{
				split_obtuse_boundary(i);
			}
		}
		// reset counter (only here, not after the second loop)
		n0 = n;

		// first loop, refine all obtuse triangles
		for (int i=0; i<n0; i++)
		{
			if (is_obtuse(i,abstol))
			{
				//System.out.println("Splitting obutse element " + i);
				split_red_yellow(i);
			} // if
		} // for i
		// second loop, refine all remainig triangles
		for (int i=0; i<n0; i++)
		{
			// internal logic avoids resplitting of already split elements
			split(i);
		} // for i
	} // refine_nonobtuse

	public final void split_obtuse_boundary(int t1_)
	{
		// determine the longest side
		// check if this is a boundary, if yest, split by falling the plumb line

		final int t1 = t1_+1;

		// search for the neighbour opposit the obtuse angle
		double lmax  = 0;
		int    imax  = 0;

		for (int i = 0; i<3; i++)
		{
			int    p1 = get_vertex(t1-1,(i+1) % 3);
			double x1 = E_tree.P[p1-1][0];
			double y1 = E_tree.P[p1-1][1];
			int    p2 = get_vertex(t1-1,(i+2) % 3);
			double x2 = E_tree.P[p2-1][0];
			double y2 = E_tree.P[p2-1][1];
			double l  = Math.hypot(x2-x1,y2-y1);
			if (l > lmax)
			{
				imax = i;
				lmax = l;
			}		
		}
		int t2 = get_neighbour(t1-1,imax);

	        if (t2 < 0)
		{
			int p0 = get_vertex(t1-1, imax);
			int p1 = get_vertex(t1-1, (imax+1) % 3);
			int p2 = get_vertex(t1-1, (imax+2) % 3);
			double x0 = E_tree.P[p0-1][0];
			double y0 = E_tree.P[p0-1][1];
			double x1 = E_tree.P[p1-1][0];
			double y1 = E_tree.P[p1-1][1];
			double x2 = E_tree.P[p2-1][0];
			double y2 = E_tree.P[p2-1][1];
			
			// point on boundary
			double [] xy3 = plumb_line(x1,y1,x2,y2,x0,y0);
			int p3 = E_tree.add_point(xy3[0],xy3[1]);
			// add both new triangles
			int c1 = add_child(p1,p0,p3,0,0,0,t1);
			set_child(t1-1,0,c1);
			int c2 = add_child(p2,p0,p3,0,0,0,t1);
			set_child(t1-1,1,c2);			
			// TODO split boundary
		} // if t2 < 0
	} // split_obtuse_boundary

	// reference Korotov 2002
	// TODO this requires only red and yellow triangles in the input mesh
	private final void split_red_yellow(final int t1_)
	{
// TODO only split if not yet split
		final int t1 = t1_+1;

		// search for the neighbour opposit the obtuse angle
		double lmax  = 0;
		int    imax  = 0;
//		int t2 = 0;
//		p1 = get_vertex(t1-1,1);
//		p2 = get_vertex(t1-1,2);

		for (int i = 0; i<3; i++)
		{
			int    p1 = get_vertex(t1-1,(i+1) % 3);
			double x1 = E_tree.P[p1-1][0];
			double y1 = E_tree.P[p1-1][1];
			int    p2 = get_vertex(t1-1,(i+2) % 3);
			double x2 = E_tree.P[p2-1][0];
			double y2 = E_tree.P[p2-1][1];
			double l  = Math.hypot(x2-x1,y2-y1);
			if (l > lmax)
			{
				imax = i;
				lmax = l;
			}		
		}
		    int t2 = get_neighbour(t1-1,imax);

		    if (t2 < 0)
		    {
		    	// neighbour is a boundary, apply 7-split
		    	// TODO this is more complicated if the bnd triangly is strongly obtuse
		    	// split_7(t1,
		    } else {
		    // verify, that the neighbour is not yet refined	
		    if (!is_leaf(t2-1))
		    {
		    	//throw new RuntimeException("Two obtuse elements bordering triangle " + t2 + " at their longest edge");
		    	System.err.println("Two obtuse elements bordering triangle " + t2 + " at their longest edge");
			return;
		    }
		    // TODO, verify, that this neighbour is not obtuse by itself

		    // the acute neighbour is ABC
		    // the obtuse triangle is ADB

		    int pa  = get_vertex(t1-1,(imax+1) % 3);
		    int pb  = get_vertex(t1-1,(imax+2) % 3);
		    int pd  = get_vertex(t1-1,imax);
		    int edb =   get_edge(t1-1,(imax+1) % 3);
		    int ead =   get_edge(t1-1,(imax+2) % 3);

		    // determine point c and non-facing edges
		    int pc  = -1;
		    int eac = -1;
		    int ebc = -1;
		    for (int i=0; i<3; i++)
		    {
			if (pa == get_vertex(t2-1,i))
			{
				if (pb == get_vertex(t2-1,(i+1) % 3))
				{
					// cw
					pc  = get_vertex(t2-1,(i+2) % 3);
					eac = get_edge(t2-1,(i+1) % 3);
					
				} else {
					// ccw
					pc  = get_vertex(t2-1,(i+1) % 3);
					eac = get_edge(t2-1,(i+2) % 3);
				}
				ebc = get_edge(t2-1,i);
				break;
			}
		    }

		    // fetch coordinates
		    double xa = E_tree.P[pa-1][0];
		    double ya = E_tree.P[pa-1][1];
		    double xb = E_tree.P[pb-1][0];
		    double yb = E_tree.P[pb-1][1];
		    double xc = E_tree.P[pc-1][0];
		    double yc = E_tree.P[pc-1][1];
		    double xd = E_tree.P[pd-1][0];
		    double yd = E_tree.P[pd-1][1];
		     
		    double a2 = xa*xa + ya*ya;
		    double b2 = xb*xb + yb*yb;
		    double c2 = xc*xc + yc*yc;
		    double d2 = xd*xd + yd*yd;

		    // coordinates of the centre-point
		    double xs = 0.5*((yb-yd)*(a2-d2) + (yd-ya)*(b2-d2))
                                   /((xa-xd)*(yb-yd) - (xb-xd)*(ya-yd));
		    double ys = 0.5*((xd-xb)*(a2-d2) + (xa-xd)*(b2-d2))
                                   /((xa-xd)*(yb-yd) - (xb-xd )*(ya-yd));
		    // TODO, test that xs and ys are convex in abc

		    // add the centre point
		    int ps = E_tree.add_point(xs,ys);
		    // TODO add interior edges
		    // not necessary?, does Tree_1d needs to store edges at all?

		    // split the side edges, but not the facing edge
		    //int e1 = get_edge(i,0);
		    int pk = E_tree.split(eac-1);
		    //int e2 = get_edge(i,1);
		    int pl = E_tree.split(ebc-1);
		    //int e3 = get_edge(i,2);
		    int pp = E_tree.split(ead-1);
		    //
		    int pq = E_tree.split(edb-1);

		    // children of first triangle
		    // TODO add_child should invoke set_child automatically
		    int c1 = add_child(pa,ps,pk,0,0,0,t1);
		    set_child(t1-1, 0, c1);
		    int c2_ = add_child(ps,pb,pl,0,0,0,t1);
		    set_child(t1-1, 1, c2_);
		    int c3 = add_child(pl,pc,pk,0,0,0,t1);
		    set_child(t1-1, 2, c3);
		    int c4 = add_child(pk,ps,pl,0,0,0,t1);
		    set_child(t1-1, 3, c4);
		    // children of second triangle
		    int c5 = add_child(pa,pp,ps,0,0,0,t2);
		    set_child(t2-1, 0 , c5);
		    int c6 = add_child(pp,pd,ps,0,0,0,t2);
		    set_child(t2-1, 1 ,c6);
		    int c7 = add_child(pd,pq,ps,0,0,0,t2);
		    set_child(t2-1, 2, c7);
		    int c8 = add_child(pq,pb,ps,0,0,0,t2);
		    set_child(t2-1, 3, c8);
		} // neighbour is not a boundary
	} // split_red_yellow

	public boolean [] is_obtuse(final double abstol)
	{
		boolean [] retval = new boolean[n];
		for (int i=0; i<n; i++)
		{
			retval[i] = is_obtuse(i, abstol);
		}
		return retval;
	} //

	// expects zero based index
	public boolean is_obtuse(final int t, final double abstol)
	{
		double [] angle = angle(t);
		for (int i=0; i<3; i++)
		{
			if (angle[i] < abstol) return true;
		}
		return false;
	} // is_obtuse

	public double [][] angle()
	{
		double [][] angle = new double[n][];
		for (int i=0; i<n; i++)
		{
			angle[i] = angle(i);
		}
		return angle;
	}

	// zero based index
	public double [] angle(final int t)
	{
		double [] angle = new double[3];
		// a*b/(||a||||b||) = cos theta > 0	
		// TODO already return obtuse edge
		double [] x = new double[3];
		double [] y = new double[3];
		for (int i=0; i<3; i++)
		{
			int pi = get_vertex(t,i);
			x[i] = E_tree.P[pi-1][0];
			y[i] = E_tree.P[pi-1][1];
		}
		double [] dx = new double[3];
		double [] dy = new double[3];
		int [] a = {1,2,0};
		for (int i=0; i<3; i++)
		{
			//dx[i] = x[(i+1)%3]-x[i];
			//dy[i] = y[(i+1)%3]-y[i];
			dx[i] = x[a[i]]-x[i];
			dy[i] = y[a[i]]-y[i];
		}

		int [] b = {2,0,1};
		for (int i=0; i<3; i++)
		{
			// test
			//if (dot(dx[i],dy[i],-dx[(i+1)%3],-dy[(i+1)%3]) < -abstol)
			angle[i] = dot(dx[a[i]],dy[a[i]],-dx[b[i]],-dy[b[i]]);
		}
		return angle;
	} // angle int t

	public static final double dot(final double x1, final double y1, 
					final double x2, final double y2)
	{
			// TODO devide individually to avoid overflow
			double dot = ((x1*x2)+(y1*y2))/(Math.hypot(x1,y1)*Math.hypot(x2,y2));
			//Math.sqrt((x1*x1 + y1*y1)*(x2*x2 + y2*y2));
			return dot;
	}

	// test function
	public final void check_neighbour()
	{
		boolean err = false;
		Hashtable<Key2,Integer> H = new Hashtable<Key2,Integer>();
		// verifies neighbourhood links
		for (int i=0; i<n; i++)
		{
		//System.out.println(i+1 + " " + is_leaf(i) + " " + get_child(i,0));
		if (is_leaf(i))
		{
			Key2 key;	
			Integer val;
			// side opposit p1
			key = new Key2(	get_vertex(i,1), get_vertex(i,2));
			val = H.get(key);
			if ( null == val )
			{
				//System.out.println("adding " + (i+1) + " " + get_vertex(i,1) + " " + get_vertex(i,2));
				H.put(key,i+1);
			} else {
				if (val != get_neighbour(i,0))
				{
					System.err.println("T " + (i+1) + " 1 " + val + " " + get_neighbour(i,0));
					err = true;
				}
			}
			// side opposit p2
			key = new Key2(	get_vertex(i,0), get_vertex(i,2));
			val = H.get(key);
			if ( null == val )
			{
				//System.out.println("adding " + (i+1) + " " + get_vertex(i,0) + " " + get_vertex(i,2));
				H.put(key,i+1);
			} else {
				if (val != get_neighbour(i,1))
				{
					System.err.println("T " + (i+1) + " 2 " + val + " " + get_neighbour(i,1));
					err = true;
				}
			}
			// side opposit p3
			key = new Key2(	get_vertex(i,0), get_vertex(i,1));
			val = H.get(key);
			if ( null == val )
			{
				//System.out.println("adding " + (i+1) + " " + get_vertex(i,1) + " " + get_vertex(i,2));
				H.put(key,i+1);
			} else {
				if (val != get_neighbour(i,2))
				{
					//System.out.println("failing " + (i+1) + " " + get_vertex(i,1) + " " + get_vertex(i,2));
					System.err.println("T " + (i+1) + " 3 " + val + " " + get_neighbour(i,2));
					err = true;
				}
			}
		} // is_leave
		} // for i
		// check the boundary
		// boundaries are not explictely saved, but if any neigh
		//for (int

		if (err) throw new RuntimeException("Neighbourhood is not consistent");
	} // check_neighbour()

	// 0 : not contained
	// 1..n : contained in triangle i
	public final int [] isin2(final double [][] x0)
	{
		int [] in = new int[x0.length];

		for (int i=0; i<x0.length; i++)
		{
			in[i] = isin(x0[i]);
		}
		System.out.println(x0.length);
		System.out.println(in.length);
		return in;
	} // isin double [][]

	// 0 : not contained
	// 1..n : contained in triangle i
	public final int isin(final double [] x0)
	{
		// recursively search all roots (triangles in the initial mesh)
		int i=0;
		int retval = 0;
		while (    (i < n)              // still more triangles available
                        && (0 == retval)        // containing triangle not yet found
                        && (0 == get_parent(i)) // this is still a root
                      )
		{
			retval = isin(i, x0);
			i++;
		}
		return retval;
	} // isin double []

	private final int isin(final int i, final double [] x0)
	{
		// check if triangle contains the point
		int retval = 0;
		//System.out.println(i);
		if ( contains(i, x0) )
		{
			//System.out.println(i + " " + is_leaf(i));
			if (is_leaf(i))
			{
				// set return value to one based index of current triangle
				retval = i+1;
			} else {
			for (int j = 0; j < 4; j++)
			{
				final int child = get_child(i, j);
				retval    = isin(child-1,x0);
				if (retval > 0)
				{
					break;
				} // if reval > 0
			} // for j
			} // else
		} // if retval > 0
		return retval;
	} // isin in double []

	private final boolean contains(final int i, final double [] x0)
	{
		// fetch coordinates
		final double [][] P = E_tree.get_P();
		final double [][] x = new double [3][2];
		for (int j=0; j<3; j++)
		{
			for (int k=0; k<2; k++)
			{
				x[j][k] = P[get_vertex(i,j)-1][k];
			}
		}
//		double x1 = P[get_vertex(i,0)-1][0];
//		double y1 = P[get_vertex(i,0)-1][1];
//		double x2 = P[get_vertex(i,1)-1][0];
//		double y2 = P[get_vertex(i,1)-1][1];
//		double x3 = P[get_vertex(i,2)-1][0];
//		double y3 = P[get_vertex(i,2)-1][1];
		return contains(x,x0);
	} // contains

	public static final boolean contains(final double [][] x, final double [] x0)
	{
		// regression vector
		double b1 = x0[0] - x[2][0];
		double b2 = x0[1] - x[2][1];
		// regression matrix
		double a11 = x[0][0] - x[2][0]; double a12 = x[1][0]-x[2][0];
		double a21 = x[0][1] - x[2][1]; double a22 = x[1][1]-x[2][1];
		// invert 2x2 matrix
		double idet = 1.0/(a11*a22 - a12*a21);
		double ia11 =  idet*a22; double ia12 = -idet*a12;
		double ia21 = -idet*a21; double ia22 =  idet*a11;
		// solve system for barycentric coordinates
		double p1 = ia11*b1 + ia12*b2;
		double p2 = ia21*b1 + ia22*b2;
		double p3 = 1.0 - p1 - p2;
//		System.out.println(p1 + " " + p2);
		// check that coordinates are convex
		boolean c = (p1 >= 0) && (p2 >= 0) && (p3 >= 0)
		         && (p1 <= 1) && (p2 <= 1) && (p3 <= 1);	
		return c;
	}

	public void fill_R()
	{
		R = new int[n];
		for (int i=0; i<n; i++)
		{
			R[i] = i;
		}
	}

	public static void ArrayPrint2d(final int [][] a)
	{
		for (int i=0; i<a.length; i++)
		{
			System.out.println(Arrays.toString(a[i]));
		}
	}


	public static final double [] plumb_line(final double x1, final double y1,
			final double x2, final double y2, final double x0, final double y0)
	{
		double p =   ((x0-x2)*(x1-x2) + (y0 - y2)*(y1 - y2))
	    			/ ( (x1-x2)*(x1-x2) + (y1 - y2)*(y1-y2) );
		double [] ret = new double[2];
		ret[0] = p*x1 + (1-p)*x2;
		ret[1] = p*y1 + (1-p)*y2;
		return ret;
	} // plumb_line
} // class Tree_2d

