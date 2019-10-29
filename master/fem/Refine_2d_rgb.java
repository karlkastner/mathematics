// Wed Jun 20 23:01:14 MSK 2012
// Karl Kästner, Berlin

import java.util.Hashtable;

class Refine_2d_rgb
{
		// member variables
		private Mesh_2d mesh;
		private Hashtable<Key2,int[]> Sh;
		private Hashtable<Key2,Integer> Ph;
		private double L[][];
		private int Lx[];
		private int af[];
		final double EPS=1e-16;


	// constructor
	public void Refine_2d_rgb(Mesh_2d cmesh)
	{
		mesh = cmesh;
		Sh   = new Hashtable<Key2,int[]>(mesh.nt);
		Ph   = new Hashtable<Key2,Integer>(mesh.np);
		af   = new int[mesh.nt];
		L    = new double[mesh.nt][3];
		Lx   = new int[mesh.nt];
	} // Refine_2d_rgb

	public void refine(final int M[])
	{
		for (int idx=0; idx< mesh.nt; idx++)
		{
			// compute side length'
			L[idx][0] = Math.sqrt( (mesh.P[mesh.T[idx][1]][0] - mesh.P[mesh.T[idx][2]][0])*(mesh.P[mesh.T[idx][1]][0] - mesh.P[mesh.T[idx][2]][0])
				        + (mesh.P[mesh.T[idx][1]][1] - mesh.P[mesh.T[idx][2]][1])*(mesh.P[mesh.T[idx][1]][1] - mesh.P[mesh.T[idx][2]][1]) );
			L[idx][1] = Math.sqrt( (mesh.P[mesh.T[idx][0]][0] - mesh.P[mesh.T[idx][2]][0])*(mesh.P[mesh.T[idx][0]][0] - mesh.P[mesh.T[idx][2]][0])
				        + (mesh.P[mesh.T[idx][0]][1] - mesh.P[mesh.T[idx][2]][1])*(mesh.P[mesh.T[idx][0]][1] - mesh.P[mesh.T[idx][2]][1]) );
			L[idx][2] = Math.sqrt( (mesh.P[mesh.T[idx][0]][0] - mesh.P[mesh.T[idx][1]][0])*(mesh.P[mesh.T[idx][0]][0] - mesh.P[mesh.T[idx][1]][0])
				        + (mesh.P[mesh.T[idx][0]][1] - mesh.P[mesh.T[idx][1]][1])*(mesh.P[mesh.T[idx][0]][1] - mesh.P[mesh.T[idx][1]][1]) );

			// find longest sides
			// TODO: this leads to non optimal results if the
			// longest side is not uinque
			int max = 0;
			if (L[idx][1] > L[idx][0]) max = 1;
			if (L[idx][2] > L[idx][max]) max = 2;
			Lx[idx] = max;
		} // for i

		// push boundaries
		for (int bdx=0; bdx<mesh.nb; bdx++)
		{
			Key2 key2 = new Key2(mesh.Bc[bdx][0], mesh.Bc[bdx][1]);
			int [] val = new int[4];
			val[0] = -(bdx+1);
			val[1] = 0;
			Sh.put(key2, val);
		} // for bdx
		
	
		// push triangle sides
		// TODO: push sides of old triangles here
		for (int tdx=0; tdx<mesh.nt; tdx++)
		{
			Key2 key2 = new Key2(mesh.T[tdx-1][0], mesh.T[tdx-1][1]);
			int [] val = new int[4];
			val[0] = tdx; val[1] = 3;
			val = Sh.put(key2,val);
			if (null != val)
			{
				val[2] = tdx; val[3] = 3;
				Sh.put(key2, val);
			}

			key2 = new Key2(mesh.T[tdx-1][0], mesh.T[tdx-1][2]);
			val = new int[4];
			val[0] = tdx; val[1] = 2;
			val = Sh.put(key2,val);
			if (null != val)
			{
				val[2] = tdx; val[3] = 2;
				Sh.put(key2, val);
			}

			key2 = new Key2(mesh.T[tdx-1][1], mesh.T[tdx-1][2]);
			val = new int[4];
			val[0] = tdx; val[1] = 1;
			val = Sh.put(key2,val);
			if (null != val)
			{
				val[2] = tdx; val[3] = 1;
				Sh.put(key2,val);
			}
		} // for tdx

		
		// split edges for each triangle to be refined
		for (int mdx=0; mdx<M.length; mdx++)
		{
			split_edges(M[mdx]);
		} // for mdx

/*
//		keySet = Sh.keySet();
//		K = keySet.toArray();
//		for kdx = 1:length(K)
//			key = K(kdx);
//			% remove hanging node from hash
//			val = Sh.get(key)'
//		}
*/

		// for each triangle
		// TODO can be simplified by noting the affected elements
		for (int fdx=0; fdx<af.length; fdx++)
		{
			int tdx=af[fdx];
			if (tdx > 0) partition(tdx);
		} // for fdx

	} // refine_2d

	private void split_edges(final int tdx)
	{
		af[tdx-1] = 1;
		int p1 = mesh.T[tdx-1][0];
		int p2 = mesh.T[tdx-1][1];
		int p3 = mesh.T[tdx-1][2];

		// remove edges
		remove_edges(tdx);

		// split its three sides
		Key2 key2 = new Key2(p2,p3);
		Integer val = Ph.get(key2);
		if (null != val)
		{
			// TODO : midpoint function
			// TODO : addpoint/getpoint function
			mesh.P[mesh.np][0] = 0.5*(mesh.P[p2][0] + mesh.P[p3][0]);
			mesh.P[mesh.np][1] = 0.5*(mesh.P[p2][1] + mesh.P[p3][1]);
			mesh.np = mesh.np+1;

			Ph.put(key2,mesh.np);
			// recursively split the longest side of the neighbouring triangle
			if (mesh.N[tdx-1][0] > 0)
			{
				split_longest(mesh.N[tdx-1][0]);
			} else {
				split_boundary(-mesh.N[tdx-1][0], mesh.np);
			}
		}

		key2 = new Key2(p1,p3);
		val = Ph.get(key2);
		if (null != val)
		{
			mesh.P[mesh.np][0] = 0.5*(mesh.P[p1][0] + mesh.P[p3][0]);
			mesh.P[mesh.np][1] = 0.5*(mesh.P[p1][1] + mesh.P[p3][1]);
			mesh.np = mesh.np+1;
			Ph.put(key2,mesh.np);
			if (mesh.N[tdx-1][1] > 0)
			{
				split_longest(mesh.N[tdx-1][1]);
			} else {
				split_boundary(-mesh.N[tdx-1][1], mesh.np);
			}
		}

		key2 = new Key2(p1,p2);
		val = Ph.get(key2);
		if (null != val)
		{
			mesh.np = mesh.np+1;
			mesh.P[mesh.np][0] = 0.5*(mesh.P[p1][0] + mesh.P[p2][0]);
			mesh.P[mesh.np][1] = 0.5*(mesh.P[p1][1] + mesh.P[p2][1]);
			Ph.put(key2,mesh.np);
			if (mesh.N[tdx-1][2] > 0)
			{
				split_longest(mesh.N[tdx-1][2]);
			} else {
				split_boundary(-mesh.N[tdx-1][2], mesh.np);
			}
		}
	} // split_edges()
	
	private void split_longest(final int tdx)
	{
		// mark triangle as processed
		af[tdx-1] = 1;

		int p1, p2, n3;

		// remove the edges
		remove_edges(tdx);

		// get longest side
		switch(Lx[tdx])
		{
			case 0 :
				p1 = mesh.T[tdx-1][1];
				p2 = mesh.T[tdx-1][2];
				n3 = mesh.N[tdx-1][1];
			break;
			case 1 :
				p1 = mesh.T[tdx-1][2];
				p2 = mesh.T[tdx-1][0];
				n3 = mesh.N[tdx-1][2];
			break;
			case 2 :
				p1 = mesh.T[tdx-1][0];
				p2 = mesh.T[tdx-1][1];
				n3 = mesh.N[tdx-1][3];
			break;
			default :
				// should not happen
				throw(new RuntimeException());
		} // switch

		Key2 key2 = new Key2(p1,p2);
		Integer val = Ph.get(key2);
		// longest side was not yet split
		if (null != val)
		{
			mesh.np = mesh.np+1;
			mesh.P[mesh.np][0] = 0.5*(mesh.P[p1][0] + mesh.P[p2][0]);
			mesh.P[mesh.np][1] = 0.5*(mesh.P[p1][1] + mesh.P[p2][1]);
			Ph.put(key2,mesh.np);
			if (n3 > 0)
			{
				// recursive splitting of longest edge
				split_longest(n3);
			} else {
				split_boundary(-n3, mesh.np);
			}
		}
	} // split_longest()
	
	private void split_boundary(final int bdx, final int pc)
	{
		mesh.Bc[mesh.nb][0] = pc;
		mesh.Bc[mesh.nb][1] = mesh.Bc[bdx][1];
	//	Bc(mesh.nb,3) = Bc(bdx,3);
		mesh.Bc[bdx-1][1] = pc;
		// push sides
		int [] val = new int[4];
		val[0] = -bdx;
		Sh.put(new Key2(mesh.Bc[bdx-1][0],mesh.Bc[bdx-1][1]),val);
		val = new int[4];
		val[0] = -(mesh.nb+1);
		Sh.put(new Key2(mesh.Bc[mesh.nb][0],mesh.Bc[mesh.nb][1]),val);
		mesh.nb = mesh.nb+1;
	} // split_boundary()
	
	private void remove_edges(final int tdx)
	{
		// remove old edge references
		// TODO, this has to go into split_edge and split_longest not into partition
		// avoid repitative removal of same entry
		// reason: old reference might be used, if one triangle is still ond and the other is new
		Key2 key2 = new Key2(mesh.T[tdx-1][0], mesh.T[tdx-1][1]);
		int [] val = Sh.remove(key2);
		// TODO change here
		//if (length(val) > 2)
		if (val[2] > 0)
		{
			if (tdx == val[0])
			{
				val[0] = val[2]; val[1] = val[3];
			}
			val[2] = 0; val[3] = 0;
			Sh.put(key2,val);
		} else {
 			// already removed
			/// TODO problem here
			if (val[0] > 0 && tdx != val[0])
			{
				Sh.put(key2,val);
			}
		}

		key2 = new Key2(mesh.T[tdx-1][0], mesh.T[tdx-1][2]);
		val = Sh.remove(key2);
		if (val[2] > 0)
		{
			if (tdx == val[0])
			{
				val[0] = val[2]; val[1] = val[3];
			}
			val[2] = 0; val[3] = 0;
			Sh.put(key2,val);
		} else {
			// TODO
			if (val[0] > 0 && tdx != val[0])
			{
				Sh.put(key2,val);
			}
		}

		key2 = new Key2(mesh.T[tdx-1][1], mesh.T[tdx-1][2]);
		val = Sh.remove(key2);
		if (val[2] > 0)
		{
			if (tdx == val[0])
			{
				val[0] = val[2]; val[1] = val[3];
			}
			val[2] = 0; val[3] = 0;
			Sh.put(key2,val);
		} else {
			if (val[0] > 0 && tdx != val[0])
			{
				Sh.put(key2,val);
			}
		}
	} // remove_edges()
	
	private void partition(final int tdx)
	{	
		int p1, p2, p3;
		double l1,l2,l3;
		switch (Lx[tdx-1])
		{
			case 1 :
				p1 = mesh.T[tdx-1][0];
				p2 = mesh.T[tdx-1][1];
				p3 = mesh.T[tdx-1][2];
				l1=L[tdx][0];
				l2=L[tdx][1];
				l3=L[tdx][2];
				break;
			case 2 :
				p1 = mesh.T[tdx-1][1];
				p2 = mesh.T[tdx-1][2];
				p3 = mesh.T[tdx-1][0];
				l1=L[tdx][1];
				l2=L[tdx][2];
				l3=L[tdx][0];
				break;
			case 3 :
				p1 = mesh.T[tdx-1][2];
				p2 = mesh.T[tdx-1][0];
				p3 = mesh.T[tdx-1][1];
				l1=L[tdx][2];
				l2=L[tdx][1];
				l3=L[tdx][0];
				break;
			default :
				// should never happen
				throw(new RuntimeException());
		} // swich
	
		// check, wether this triangle is obtuse
		double cos_alpha = (l2*l2 + l3*l3 - l1*l1)/(2*l2*l3);
		boolean obtuse = (cos_alpha+Math.sqrt(EPS) < 0);
	
		// fetch new points
		// TODO, use the T(4:6) and also preallocate
		//p12 = T(idx,6);
		//p13 = T(idx,5);
		//p23 = T(idx,4);
		Integer p12 = Ph.get(new Key2(p1,p2));
		Integer p13 = Ph.get(new Key2(p1,p3));
		Integer p23 = Ph.get(new Key2(p2,p3));
	
		// determine, which sides in addition to 1-2 are split
		int s = (null != p23 ? 1:0) + 2*(null != p13 ? 1:0) + 4*(null != p12 ? 1:0);
		// partition the triangle according to 4 different cases
		switch (s)
		{
			case 0 :
				// no side split
				break;
			case 1 : // only longest side was split
				add_triangle(tdx,p1,p2,p23);
				mesh.nt = mesh.nt+1;
				add_triangle(mesh.nt,p1,p23,p3);
				break;
			case 3 : // longest and right side were split
				add_triangle(tdx,p1,p2,p23);
				mesh.nt = mesh.nt+1;
				add_triangle(mesh.nt,p1,p23,p13);
				mesh.nt = mesh.nt+1;
				add_triangle(mesh.nt,p23,p3,p13);
				break;
			case 6 : // longest and left side were split
				add_triangle(tdx,p1,p12,p23);
				mesh.nt = mesh.nt+1;
				add_triangle(mesh.nt,p2,p23,p12);
				mesh.nt = mesh.nt+1;
				add_triangle(mesh.nt,p1,p23,p3);
				break;
			case 7 : // regular refinement, longest and both other sides were split
			if (!obtuse)
			{
				add_triangle(tdx,p12,p13,p23);
				mesh.nt = mesh.nt+1;
				add_triangle(mesh.nt,p1,p12,p13);
				mesh.nt = mesh.nt+1;
				add_triangle(mesh.nt,p2,p12,p23);
				mesh.nt = mesh.nt+1;
				add_triangle(mesh.nt,p3,p13,p23);
			} else {
				add_triangle(tdx,p1,p12,p23);
				mesh.nt = mesh.nt+1;
				add_triangle(mesh.nt,p2,p23,p12);
				mesh.nt = mesh.nt+1;
				add_triangle(mesh.nt,p1,p23,p13);
				mesh.nt = mesh.nt+1;
				add_triangle(mesh.nt,p23,p3,p13);
			}
			break;
			default:
				throw(new RuntimeException("longest side was not split"));
		} // switch
	} // function partition
	
	private void add_triangle(final int tdx, final int p1, final int p2, final int p3)
	{
		// TODO, add triangle should be part of Mesh_2d
		mesh.T[tdx-1][0] = p1;
		mesh.T[tdx-1][1] = p2;
		mesh.T[tdx-1][2] = p3;

		// neighbour 1
		Key2 key2 = new Key2(p2,p3);
		int [] val = Sh.remove(key2);
		if (null == val)
		{
			val = new int[4];
			val[0] = tdx; val[1] = 1;
			Sh.put(key2,val);
		} else {
			if (val[0] > 0)
			{
				// neighbour is a triangle
				mesh.N[tdx-1][0] = val[0];
				mesh.N[val[0]-1][val[2]] = tdx;
			} else {
				// neighbour is a boundary
				mesh.N[tdx-1][0] = val[0];
			}
		}
		// neighbour 2
		key2 = new Key2(p1,p3);
		val = Sh.remove(key2);
		if (null == val)
		{
			val = new int[4];
			val[0] = tdx; val[1] = 2;
			Sh.put(key2,val);
		} else {
			if (val[0] > 0)
			{
				// neighbour is a triangle
				mesh.N[tdx-1][1] = val[0];
				mesh.N[val[0]-1][val[1]] = tdx;
			} else {
				// neighbour is a boundary
				mesh.N[tdx-1][1] = val[0];
			}
		}
		// neighbour 3
		key2 = new Key2(p1,p2);
		val = Sh.remove(key2);
		if (null == val)
		{
			val = new int[4];
			val[0] = tdx; val[1] = 3;
			Sh.put(key2, val);
		} else {
			if (val[0] > 0)
			{
				// neighbour is a triangle
				mesh.N[tdx-1][2] = val[0];
				mesh.N[val[0]-1][val[2]] = tdx;
			} else {
				// neighbour is a boundary
				mesh.N[tdx-1][2] = val[0];
			}
		}
	} // add_triangle
	
} // class Refine_2d_rgb

