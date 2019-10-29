// Tue Jul 24 23:49:50 MSK 2012
// Karl Kästner, Berlin

// todo, multigrid setup, parallel and distributed functions
import Jama.*;
import java.util.Hashtable;

// make common interface for meshes
public final class Mesh_3d
{
	// point coordinates
	// gets extended by refinement
	public double [][] P;
	// number of points
	public int np;
	// conforming "red" and "green" tetrahedra array generated from tree
	// columns get extended by promotion
	public int [][] T;
	// alternative triagulation with basis function p/2 for error estimation
	public int [][] S;
	public int ns;
	// number of triangles
	public int nt;
	// boundary faces
	public int [][] Bc;
	// element neighbours
	public int [][] N;
	// number of boundary faces
	public int nb;
	// dimension
	public static final int DIM = 3;
	// centre coordinates
	public double [][]	C;
	// maximum side length
	public double []	h_side;
	// degeneracy of the element (volume / h_side^3)
	public double []	degen;
	// determinant of the element
	public double []	determinant;
	// threshold for refinement
	final static double THRESH = 0.5;	

	private double P_local[][][];
	private double Phi[][][];

	// constructor
	public Mesh_3d(final int np, final int nt, final int nb)
	{
		P = new double[np][DIM];
		this.np = 0;
		T = new int[nt][DIM+1];
		this.nt = 0;
		Bc = new int[nb][DIM];
		this.nb = 0;
	} // constructor

	public Mesh_3d(double [][] P, int [][] T, int [][] Bc)
	{
		this.P = P;
		this.np = P.length;
		this.T = T;
		this.nt = T.length;
		this.Bc = Bc;
		this.nb = Bc.length;
		//this.Phi = null;
	}

	public final void add_T(final int p0, final int p1, final int p2, final int p3)
	{
		T[nt][0] = p0;
		T[nt][1] = p1;
		T[nt][2] = p2;
		T[nt][3] = p3;
		nt++;
	} // add_T()

	public final void add_Bc(final int p0, final int p1, final int p2)
	{
		Bc[nb][0] = p0;
		Bc[nb][1] = p1;
		Bc[nb][2] = p2;
		if (p0 == p1 || p0 == p2 || p1 == p2) throw new RuntimeException();
		//if (P(p0 == p1 || p0 == p2 || p1 == p2) throw new RuntimeException();
		nb++;
	} // add_Bc()

	// determine degree of basis function polynomial
	//  1,4,10,20,35,... -> 0,1,2,3,4,...
	// lt2 = 1/6*nv^3 + nv^2 + 11/6*nv + 1
	private final int get_nv()
	{
		int nt2 = T[0].length; // TODO, dangerous if T empty
		double q = Math.pow(3.0*nt2+Math.sqrt(9.0*nt2*nt2-1.0/27.0),1.0/3.0);
		// the result is only an integer, round only used
		// to not to cast accidently down a number perturbed by finite precission
 		int nv = (int) Math.round(1.0/(3.0*q) + q - 2.0);
		return nv;
	} // get nv
	
/*
	public double [][] P() { return T_tree.P; }
	public int [][] T() { return T_tree.T; }
	public int [][] Bc() { return T_tree.Bc; }
	public int np() { return T_tree.np; }
	public int nt() { return T_tree.nt; }
	public int nb() { return T_tree.nb; }
*/

	// what is faster, no prefetch? prefetch only point coordinates? prefetch test functions?
	// locally promote ?
	public final void prefetch()
	{
		int nt1 = nt;
		int nt2 = T[0].length; // dangerous if T empty
		// determine degree of basis function polynomial
		int nv = get_nv();

		// allocate memory
		Phi         = new double[nt1][nt2][nt2];
		// TODO, this does not really exist, or 3*4 side angles
		// s_angle     = new double[nt1][3];
		// todo, only the maximum side length is actually required
		// h_side      = new double[nt1][DIM*(DIM-1)]; // 6 in 3D
		h_side	    = new double[nt1];
		// centre coordinates
		C           = new double[nt1][DIM];
		// TODO, replace determinant by hypervolume
		determinant = new double[nt1];
		// measure of the degeneracy of the element
		degen	    = new double[nt1];
		// int 3D, this is an area
		// l_boundary  = new double[Bc.length];
		P_local     = new double[nt1][DIM+1][DIM];
		
		// temporary arrays
		double [][]   A    = new double[DIM+1][DIM+1];
		double [][]   A_   = new double[nt2][DIM];
		double [][]   Va   = new double[nt2][(nv+1)*(nv+2)*(nv+3)/6];
		Matrix        mVa  = new Matrix(Va);
		Matrix        mA   = new Matrix(A);
		
		for (int tdx=0; tdx<nt1; tdx++)
		{
			int [] T_tdx = T[tdx];

			// prefetched point coordinates (TODO, use below and do not rely on compiler)
			P_local[tdx][0][0] = P[T_tdx[0]-1][0]; P_local[tdx][0][1] = P[T_tdx[0]-1][1]; P_local[tdx][0][2] = P[T_tdx[0]-1][2];
			P_local[tdx][1][0] = P[T_tdx[1]-1][0]; P_local[tdx][1][1] = P[T_tdx[1]-1][1]; P_local[tdx][1][2] = P[T_tdx[1]-1][2];
			P_local[tdx][2][0] = P[T_tdx[2]-1][0]; P_local[tdx][2][1] = P[T_tdx[2]-1][1]; P_local[tdx][2][2] = P[T_tdx[2]-1][2];
			P_local[tdx][3][0] = P[T_tdx[3]-1][0]; P_local[tdx][3][1] = P[T_tdx[3]-1][1]; P_local[tdx][3][2] = P[T_tdx[3]-1][2];

			//  prefetch triangle points
			for (int j=0; j<nt2; j++)
			{
				A_[j][0] = P[T_tdx[j]-1][0];
				A_[j][1] = P[T_tdx[j]-1][1];
				A_[j][2] = P[T_tdx[j]-1][2];
			}
			
			
			// construct the Vandermonde matrix
			FEM.vander_3d(Va, A_, nv);
			
			// invert the Vandermonde matrix
			// TODO fetch singular matrices and use SVD to compute the pseudoinverse
			Matrix mC = mVa.inverse();

			// write to output array
			Phi[tdx] = mC.getArray();

			// calculate element centre coordinate
			C[tdx][0] = 0.25*(P[T_tdx[0]-1][0] + P[T_tdx[1]-1][0] + P[T_tdx[2]-1][0] + P[T_tdx[3]-1][0] );
			C[tdx][1] = 0.25*(P[T_tdx[0]-1][1] + P[T_tdx[1]-1][1] + P[T_tdx[2]-1][1] + P[T_tdx[3]-1][1] );
			C[tdx][2] = 0.25*(P[T_tdx[0]-1][2] + P[T_tdx[1]-1][2] + P[T_tdx[2]-1][2] + P[T_tdx[3]-1][2] );

			// fetch corner points
			// TODO do not fecht again, use submatrix
			A[0][0] = 1; A[0][1] = P[T_tdx[0]-1][0]; A[0][2] = P[T_tdx[0]-1][1]; A[0][3] = P[T_tdx[0]-1][2];
	                A[1][0] = 1; A[1][1] = P[T_tdx[1]-1][0]; A[1][2] = P[T_tdx[1]-1][1]; A[1][3] = P[T_tdx[1]-1][2];
	                A[2][0] = 1; A[2][1] = P[T_tdx[2]-1][0]; A[2][2] = P[T_tdx[2]-1][1]; A[2][3] = P[T_tdx[2]-1][2];
	                A[3][0] = 1; A[3][1] = P[T_tdx[3]-1][0]; A[3][2] = P[T_tdx[3]-1][1]; A[3][3] = P[T_tdx[3]-1][2];

			// calculate element determinant
			determinant[tdx] = Math.abs(mA.det());

			// TODO: warn for clockwise (area<0) and degenerated triangles (area/h_side < eps)

			// calculate the side points
			// todo, actually only the maximum side is required, not every side length
			int k = 0;
			double [] c = new double [3];
			h_side[tdx] = 0; // paranoid, java should initialise with 0
			// TODO, no magic numbers
			// alternatives are geometric mean, 1-norm, 2-norm
			double h;
			switch (1)
			{
				case 0 :
					// maximum edge (inf norm)
					h =             FEM.dist3(P[T_tdx[0]-1], P[T_tdx[1]-1]);
					h = Math.max(h, FEM.dist3(P[T_tdx[0]-1], P[T_tdx[2]-1]));
					h = Math.max(h, FEM.dist3(P[T_tdx[0]-1], P[T_tdx[3]-1]));
					h = Math.max(h, FEM.dist3(P[T_tdx[1]-1], P[T_tdx[2]-1]));
					h = Math.max(h, FEM.dist3(P[T_tdx[1]-1], P[T_tdx[3]-1]));
					h = Math.max(h, FEM.dist3(P[T_tdx[2]-1], P[T_tdx[3]-1]));
					h_side[tdx] = h;
				break;
				case 1 :
					// diameter of circumcircle
					c[0] = 0.25*( P[T_tdx[0]-1][0] + P[T_tdx[1]-1][0] + P[T_tdx[2]-1][0] + P[T_tdx[3]-1][0] );
					c[1] = 0.25*( P[T_tdx[0]-1][1] + P[T_tdx[1]-1][1] + P[T_tdx[2]-1][1] + P[T_tdx[3]-1][1] );
					c[2] = 0.25*( P[T_tdx[0]-1][2] + P[T_tdx[1]-1][2] + P[T_tdx[2]-1][2] + P[T_tdx[3]-1][2] );
					h_side[tdx] = 2*FEM.dist3(c, P[T_tdx[0]-1]);
				break;
			} // switch (0)

			// in 3D using the angles of the faces is expensive (12 angles)
			// so the degeneracy factor volume / r_max^3 is used (must be bound away from zero for consisten meshes
			degen[tdx] = Math.pow(determinant[tdx] / (h_side[tdx]*h_side[tdx]*h_side[tdx]), 1.0/3.0);

			// sine of interior angles of faces
			// TODO
		} // for tdx (each tetra)

		// TODO separate boundary area (this is only used for test purposes)
		/*
		for (int idx=0; idx<Bc.length; idx++)
		{
			l_boundary[idx] = Math.sqrt( (P[Bc[idx][0]-1][0] - P[Bc[idx][1]-1][0])*(P[Bc[idx][0]-1][0] - P[Bc[idx][1]-1][0])
						   + (P[Bc[idx][0]-1][1] - P[Bc[idx][1]-1][1])*(P[Bc[idx][0]-1][1] - P[Bc[idx][1]-1][1]) );
		}
		*/
	} // prefetch_2d

	// compute neighbouring elements
	// TODO : generate neighbours implicitely and much faster from the tree
	public final void element_neighbours()
	{
		// allocate memory
		N = new int[nt][4];

		// hash of triangle sides
		// key: lower_point_index, higher_point_index
		// value: triangle index, index of opposit point
		// if boundary, boundary index and negative number of boundary
		Hashtable<Key3,int[]> Sh = new Hashtable<Key3,int[]>(4*T.length);

		// push the boundaries
		for (int idx=0; idx<nb; idx++)
		{
			// push row and column in N, negative, domain
			int val[] = new int[1];
			val[0] = -idx-1;
			Sh.put(new Key3(Bc[idx][0], Bc[idx][1], Bc[idx][2]), val);
		}

		final int [][] q = { 
		{ 1, 2, 3 },
		{ 0, 2, 3 },
		{ 0, 1, 3 },
		{ 0, 1, 2 } };
//		for (int jdx=0; jdx<4; jdx++)
//			System.out.println(q[jdx][0] + " " + q[jdx][1] + " " + q[jdx][2]);

		for (int idx=0; idx<nt; idx++)
		{
			for (int jdx=0; jdx<4; jdx++)
			{
				// point at opposit side
				Key3 key3 = new Key3(T[idx][q[jdx][0]], T[idx][q[jdx][1]], T[idx][q[jdx][2]]);
				int [] val = Sh.remove(key3);
				if (null == val)
				{
					val = new int[2];
					val[0] = idx;
					val[1] = jdx;
					Sh.put(key3, val);
				} else {
					// Neighbourhood relation found
					if (val[0] >= 0)
					{
						N[idx][jdx] = val[0]+1;
						N[val[0]][val[1]] = idx+1;
					} else {
						N[idx][jdx] = val[0];
					}
				}
			} // for jdx
		} // for idx
		
		if (!Sh.isEmpty())
		{
			if (DEBUG.level > 1)
			{
			Object [] K = Sh.keySet().toArray();
			System.err.println(K.length);
			for (int i=0; i<K.length; i++)
			{
				int[] val = Sh.get(K[i]);
				System.err.println(((Key3)K[i]).A[0] + " " + ((Key3)K[i]).A[1] + " " + ((Key3)K[i]).A[2]
					+ " " + val[0] + " " + N[val[0]][0] + " " + N[val[0]][1] + " " + N[val[0]][2] + " " + N[val[0]][3] );
			}
			}
			throw new RuntimeException("Mesh_3d.element_neighbours: inconsistent mesh");
		}
	} // element_neighbours()

	// assemble of mass and potential matrix
	//  quadrature coordinates and weights w,b, flag
	// TODO make the first two column of buf an integer array
	public final double [][] assemble_phi_phi( Potential_3D func, double [] w, double [][] b, int flag)
	{
		int nt1 = nt;
		int nt2 = T[0].length; // dangerous, if T is empty

		// determine degree of basis function polynomial
		int nv = get_nv();

		int mb0;
		if (0 != flag && null == func)
		{
			// diagonal mass matrix
			mb0 = nt2;
		} else {
			// non diagonal mass matrix
			// number of buffer entries per triangle
			mb0 = nt2*(nt2+1)/2;
		}

		// exact number of elements in buffer after integration
		int nb0 = nt1*mb0;

		//  allocate memory
		double [][] buf = new double[nb0][3];

		// preallocate thread local buffers
		double[][]  A   = new double[4][4];
		double[]    wfa = new double[w.length];
		double[][]  A_  = new double[nt2][3];
		double[][]  Va  = new double[nt2][(nv+1)*(nv+2)*(nv+3)/6];
		double[][]  Vq  = new double[w.length][(nv+1)*(nv+2)*(nv+3)/6];
		Matrix mA       = new Matrix(A);
		Matrix mVa      = new Matrix(Va);
		Matrix mB       = new Matrix(b);
		Matrix mVq      = new Matrix(Vq);

		Matrix mC;

		// integrate shares of each triangle
		for (int idx=0; idx<nt1; idx++)
		{
			// fetch point coordinates
			if (null == P_local)
			{
				A[0][0] = 1; A[0][1] = P[T[idx][0]-1][0]; A[0][2] = P[T[idx][0]-1][1]; A[0][3] = P[T[idx][0]-1][2];
				A[1][0] = 1; A[1][1] = P[T[idx][1]-1][0]; A[1][2] = P[T[idx][1]-1][1]; A[1][3] = P[T[idx][1]-1][2];
				A[2][0] = 1; A[2][1] = P[T[idx][2]-1][0]; A[2][2] = P[T[idx][2]-1][1]; A[2][3] = P[T[idx][2]-1][2];
				A[3][0] = 1; A[3][1] = P[T[idx][3]-1][0]; A[3][2] = P[T[idx][3]-1][1]; A[3][3] = P[T[idx][3]-1][2];
			} else {
				// values are prefetched
				A[0][0] = 1; A[0][1] = P_local[idx][0][0]; A[0][2] = P_local[idx][0][1]; A[0][3] = P_local[idx][0][2];
				A[1][0] = 1; A[1][1] = P_local[idx][1][0]; A[1][2] = P_local[idx][1][1]; A[1][3] = P_local[idx][1][2];
				A[2][0] = 1; A[2][1] = P_local[idx][2][0]; A[2][2] = P_local[idx][2][1]; A[2][3] = P_local[idx][2][2];
				A[3][0] = 1; A[3][1] = P_local[idx][3][0]; A[3][2] = P_local[idx][3][1]; A[3][3] = P_local[idx][3][2];
			}

			double volume;
			if (null == determinant)
			{
				// tetra volume
				volume = 1.0/6.0*Math.abs(mA.det());
			} else {
				volume = 1.0/6.0*determinant[idx];
			}

			// quadrature points
			Matrix mQ = mB.times(mA);

			// no prefetch
			if (null == Phi)
			{
				// fetch all coordinates, not just corner coordinates
				// todo, can be computed on the fly (except on curved boundaries)
				for (int j=0; j<nt2; j++)
				{
					A_[j][0] = P[T[idx][j]-1][0];
					A_[j][1] = P[T[idx][j]-1][1];
					A_[j][2] = P[T[idx][j]-1][2];
				}
				// vandermonde matrix
				FEM.vander_3d(Va, A_, nv);
				// get test function coefficients
				// changed from 0, 3 to 0 lt2-1
				mC = mVa.inverse();
			} else {
				// fetch precomputed test function coefficients
				mC = new Matrix(Phi[idx]);
			}

			// evaluation of the polynomial test functions at the quadrature points
			FEM.vander_3d(Vq, mQ.getMatrix(0, w.length-1, 1, 3).getArray(), nv);
			
			Matrix mPhi     = mVq.times(mC);
			double [][] phi = mPhi.getArray();

	
			// evaluate PDE-coefficient function at quadrature points
			// and premultipy values
			if (null != func)
			{
				double[] f = func.func(mQ.getMatrix(0, w.length-1, 1, 3).getArray());
	
				// premultiply values
				for (int k=0; k<w.length; k++)
				{
					wfa[k] = w[k]*f[k]*volume;
				}
			} else {
				// premultiply values
				for (int k=0; k<w.length; k++)
				{
					wfa[k] = w[k]*volume;
				}
			}


			// add share to each of the elements corner points
			int mb = mb0*idx;
			for (int adx=0; adx<nt2; adx++)
			{
				// diagonal entry
				// TODO, this can be furthermore optimised in case of diagonal mass matrices as phi_k = 1 at k=adx and zero otherwise
				double I = 0;
				for (int k=0; k<w.length; k++)
				{
					I = I + wfa[k]*phi[k][adx]*phi[k][adx];
				}
				buf[mb][0] = T[idx][adx];
				buf[mb][1] = T[idx][adx];
				buf[mb][2] = 0.5*I;
				mb++;
				if (0 == flag || null != func)
				{
					// off diagonal entries
					for (int bdx=adx+1; bdx<nt2; bdx++)
					{
						// A(p(adx),p(bdx)) = A(p(adx),p(bdx)) + I;
						I = 0;
						for (int k=0; k<w.length; k++)
						{
							I = I + wfa[k]*phi[k][adx]*phi[k][bdx];
						}
						buf[mb][0] = T[idx][adx];
						buf[mb][1] = T[idx][bdx];
						buf[mb][2] = I;
						mb++;
					} // for bdx
				} // if 0 == flag
			} // for adx
		} //  idx

		return buf;
	} // assemle_dphi_dphi()

	// assembly of stiffness matrix
	// returns buffer
	public final double [][] assemble_dphi_dphi( Potential_3D func, double [] w, double [][] b)
	{
		int nt1 = nt;
		int nt2 = T[0].length; // dangerous, if T is empty

		// get degree of basis function polynomial
 		int nv = get_nv();

		// stiffness matrix is never diagonal
		// number of buffer entries per triangle
		int mb0 = nt2*(nt2+1)/2;

		// exact number of elements in buffer after integration
		int nb0 = nt1*mb0;

		//  allocate memory
		double [][] buf = new double[nb0][3];

		// preallocate thread local buffers
		double[][]  A   = new double[DIM+1][DIM+1];
		double[]    wfa = new double[w.length];
		double[][]  A_  = new double[nt2][DIM];
		double[][]  Va  = new double[nt2][(nv+1)*(nv+2)*(nv+3)/6];
		double[][]  Vq  = new double[w.length][nv*(nv+1)*(nv+2)/6];
		double[][] dC_x = new double[nv*(nv+1)*(nv+2)/6][(nv+1)*(nv+2)*(nv+3)/6];
		double[][] dC_y = new double[nv*(nv+1)*(nv+2)/6][(nv+1)*(nv+2)*(nv+3)/6];
		double[][] dC_z = new double[nv*(nv+1)*(nv+2)/6][(nv+1)*(nv+2)*(nv+3)/6];

		// JAMA matrix wrappers for buffers
		Matrix mA       = new Matrix(A);
		Matrix mVa      = new Matrix(Va);
		Matrix mB       = new Matrix(b);
		Matrix mVq      = new Matrix(Vq);
		Matrix mDC_x    = new Matrix(dC_x);
		Matrix mDC_y    = new Matrix(dC_y);
		Matrix mDC_z    = new Matrix(dC_z);
		
		double C [][];

		// integrate shares of each triangle
		for (int idx=0; idx<nt1; idx++)
		{
			// fetch point coordinates
			if (null == P_local)
			{
				A[0][0] = 1; A[0][1] = P[T[idx][0]-1][0]; A[0][2] = P[T[idx][0]-1][1]; A[0][3] = P[T[idx][0]-1][2];
				A[1][0] = 1; A[1][1] = P[T[idx][1]-1][0]; A[1][2] = P[T[idx][1]-1][1]; A[1][3] = P[T[idx][1]-1][2];
				A[2][0] = 1; A[2][1] = P[T[idx][2]-1][0]; A[2][2] = P[T[idx][2]-1][1]; A[2][3] = P[T[idx][2]-1][2];
				A[3][0] = 1; A[3][1] = P[T[idx][3]-1][0]; A[3][2] = P[T[idx][3]-1][1]; A[3][3] = P[T[idx][3]-1][2];
			} else {
				// values are prefetched
				A[0][0] = 1; A[0][1] = P_local[idx][0][0]; A[0][2] = P_local[idx][0][1]; A[0][3] = P_local[idx][0][2];
				A[1][0] = 1; A[1][1] = P_local[idx][1][0]; A[1][2] = P_local[idx][1][1]; A[1][3] = P_local[idx][1][2];
				A[2][0] = 1; A[2][1] = P_local[idx][2][0]; A[2][2] = P_local[idx][2][1]; A[2][3] = P_local[idx][2][2];
				A[3][0] = 1; A[3][1] = P_local[idx][3][0]; A[3][2] = P_local[idx][3][1]; A[3][3] = P_local[idx][3][2];
			}

			// tetra volume
			double volume;
			if (null == determinant)
			{
				volume = 1.0/6.0*Math.abs(mA.det());
			} else {
				volume = 1.0/6.0*determinant[idx];
			}

			// quadrature points in homogeneous Cartesian coordinates
			Matrix mQ = mB.times(mA);

			// get test function coefficients
			if (null == Phi)
			{
				// fetch all point coordinates, not just for corner points
				// TODO, can be computed on the fly, except on curved boundaries
				for (int j=0; j<nt2; j++)
				{
					A_[j][0] = P[T[idx][j]-1][0];
					A_[j][1] = P[T[idx][j]-1][1];
					A_[j][2] = P[T[idx][j]-1][2];
				}
				// vandermonde matrix
				FEM.vander_3d(Va, A_, nv);

				// compute the coefficients of the polynomial test functions
				Matrix mC = mVa.inverse();
				C = mC.getArray();
			} else {
				// load precomputed test function coefficients
				C = Phi[idx];
			}

			// compute the partial derivatives of the test functions
			FEM.derivative_3d(C, dC_x, dC_y, dC_z, nv);

			// evaluate the partial derivatives of the polynomail test functions at the quadrature points
			FEM.vander_3d(Vq, mQ.getMatrix(0, w.length-1, 1, 3).getArray(), nv-1);
			Matrix mDPhi_x = mVq.times(mDC_x);
			Matrix mDPhi_y = mVq.times(mDC_y);
			Matrix mDPhi_z = mVq.times(mDC_z);
			double [][] dphi_x = mDPhi_x.getArray();
			double [][] dphi_y = mDPhi_y.getArray();
			double [][] dphi_z = mDPhi_z.getArray();

			// evaluate PDE-coefficient function at quadrature points
			// and premultipy values
			if (null != func)
			{
				double[] f = func.func(mQ.getMatrix(0, w.length-1, 1, 3).getArray());
	
				// premultiply values
				for (int k=0; k<w.length; k++)
				{
					wfa[k] = w[k]*f[k]*volume;
				}
			} else {
				// premultiply values
				for (int k=0; k<w.length; k++)
				{
					wfa[k] = w[k]*volume;
				}
			}

			// add share to each of the elements corner points
			int mb = mb0*(idx);
			for (int adx=0; adx<nt2; adx++)
			{
				// diagonal entry
				double I = 0;
				for (int k=0; k<w.length; k++)
				{
					I = I - wfa[k]*( dphi_x[k][adx]*dphi_x[k][adx]
						+ dphi_y[k][adx]*dphi_y[k][adx]   
						+ dphi_z[k][adx]*dphi_z[k][adx] );
				}
				buf[mb][0] = T[idx][adx];
				buf[mb][1] = T[idx][adx];
				buf[mb][2] = 0.5*I;
				mb++;
				// off diagonal entries
				for (int bdx=adx+1; bdx<nt2; bdx++)
				{
					// A(p(adx),p(bdx)) = A(p(adx),p(bdx)) + I;
					I = 0;
					for (int k=0; k<w.length; k++)
					{
						I = I - wfa[k]*( dphi_x[k][adx]*dphi_x[k][bdx]
							+ dphi_y[k][adx]*dphi_y[k][bdx] 
							+ dphi_z[k][adx]*dphi_z[k][bdx] );
					}
					buf[mb][0] = T[idx][adx];
					buf[mb][1] = T[idx][bdx];
					buf[mb][2] = I;
					mb++;
				} // for bdx
			} // for adx
		} //  idx

		return buf;
	} // assemble_phi_phi()


	// compute the highest nonzero derivative
	public final double [][] dV(final double [] V, final int d)
	{
		int lt1 = nt;
		int nt2 = T[0].length; // insecure
		int nv = get_nv(); // nv should equal d
		// allocate memory
		int n = (d+1)*(d+2)*(d+3)/6;

		double [][] dV      = new double[lt1][(d+1)*(d+2)/2];
		double [][] V_local = new double[n][1];
		double [][] A       = new double[nt2][4];
		double [][] Va      = new double[nt2][(nv+1)*(nv+2)*(nv+3)/6];
		Matrix      mVa     = new Matrix(Va);
		Matrix      mA      = new Matrix(A);
		Matrix mV_local     = new Matrix(V_local);

		if (V.length != np)
		{
			throw new RuntimeException("dV : length of input vector must match number of unknowns");
		}

		// compute the coefficients
		final int p = d;
        	int [] f = new int[p+1];
		f[0] = 1;
        	for (int pdx=1; pdx<=p; pdx++)
		{
                	f[pdx] = pdx*f[pdx-1];
        	}

        	final int n0 = p*(p+1)*(p+2)/6;
		int [] dC = new int[(p+1)*(p+2)/2];

		int ndx=0;
                for (int dx=p; dx>=0; dx--)
		{
			for (int dy=p-dx; dy>=0; dy--)
			{
                       		int dz=p-dx-dy;
                               	dC[ndx] = f[dx]*f[dy]*f[dz];
                               	ndx++;
			} // dy
		} // dx

		for (int idx=0; idx<lt1; idx++)
		{
			// fetch the function values
			for (int jdx=0; jdx<(d+1)*(d+2)*(d+3)/6; jdx++)
			{
				V_local[jdx][0] = V[T[idx][jdx]-1];
			}

			//  fetch triangle points
			for (int j=0; j<nt2; j++)
			{
				A[j][0] = P[T[idx][j]-1][0];
				A[j][1] = P[T[idx][j]-1][1];
				A[j][2] = P[T[idx][j]-1][2];
			}

			// construct the Vandermonde matrix
			FEM.vander_3d(Va, A, d);

			// invert the Vandermonde matrix
			// to get the test functions
			Matrix mC = mVa.inverse();

			// test/trial function coefficients
				//Matrix mC = new Matrix(Phi[idx]);
			Matrix mD = mC.times(mV_local);
			double [][] D = mD.getArray();
			
			// compute the coefficients of the partial derivatives
			for (int jdx=0; jdx<dC.length; jdx++)
			{
				dV[idx][jdx] = dC[jdx]*D[n0+jdx][0];
			}

		} // for i
		return dV;
	} // dV

	// estimates the error based on the (p+1)-th derivative
	public final Object estimate_error(final double[][] dV, final int k, final int rate)
	{
		//  get length
		int lt1 = nt;

		//  allocate memory
		double err_max = 0;
		double thresh  = 0;
		double err []  = new double[lt1];
		double nH  []  = new double[lt1];
		Object obj []  = new Object[4];

		// constants
		// C_ = [1.0/3.0 1/38 1/870 1/35430 1/882847]; % for poisson mode 0
		// C_ = [1.0/3.0, 1/13, 1/87, 1/908, 1/12360]; % for Schrödinginger mode 0
		//double [] C_ = { 1.0/9, 1.0/(3*13), 1.0/(3*87), 1.0/(3*908), 1.0/(3*12360) }; // % for Schrödinger mode 1
//		double [] C_ = { 1.0/3.0, 1.0/13.0, 1.0/87.0, 1.0/908.0, 1.0/12360.0 }; // % for Schrödinger mode 1
		double [] C_ = { 1e-0, 1e-2, 1e-4, 1e-6, 1e-8 };
//		C_[0] = 100000*C_[0];
		
		// calculate estimated norm of the second derivative per element
		for (int idx=0; idx<lt1; idx++)
		{
			// seminorm of derivativ p+1
			// for all three triangle neighbours
			// Extension to Eriksson, Johnson 1988, Adaptive Finite Element Method for Linear Elliptic Problems
			// fails for s=0 with nonsmooth data
			// no magick numbers
			for (int jdx=0; jdx<4; jdx++)
			{
				// check that neighbour exists (not a domain boundary)
				if (N[idx][jdx] > 0)
				{
					double dx      = C[idx][0] - C[N[idx][jdx]-1][0];
					double dy      = C[idx][1] - C[N[idx][jdx]-1][1];
					double dz      = C[idx][2] - C[N[idx][jdx]-1][2];
					double dr_sqr  = dx*dx + dy*dy + dz*dz;
					double nH_ = 0;
					// maximum of the partial derivatives
					for (int ddx=0; ddx<dV[0].length; ddx++)
					{
						nH_ = Math.max(nH_, Math.abs((dV[idx][ddx] - dV[N[idx][jdx]-1][ddx])*dx)/dr_sqr);
						nH_ = Math.max(nH_, Math.abs((dV[idx][ddx] - dV[N[idx][jdx]-1][ddx])*dy)/dr_sqr);
						nH_ = Math.max(nH_, Math.abs((dV[idx][ddx] - dV[N[idx][jdx]-1][ddx])*dz)/dr_sqr);
					}
					// maximum of the maximum partial derivatives of all four neighbours
					nH[idx] = Math.max(nH[idx], nH_);
				} // N(idx,jdx) > 0
			} // for jdx

			// local error of the element
			// todo, cube root of degen ?
			
			err[idx] = C_[k-2] * nH[idx] * Math.pow( h_side[idx]/degen[idx], rate);

			// global error over all elements
			err_max = Math.max(err_max, err[idx]);
		} // for idx

		// threshold for refinement
		thresh = Math.pow(THRESH, rate)*err_max;

		// stack retun values
		obj[0] = err;
		obj[1] = err_max;
		obj[2] = thresh;
		obj[3] = nH;	

		return obj;
	} // err_2d


	// Thu Jul 12 15:47:38 MSK 2012
	// Karl Kästner, Berlin
	// input p : order of the Lagrangian basis functions
	public final void promote(final int p)
	{
		if (T[0].length != 4)
		{
			throw new RuntimeException("Mesh_3d.promote: mesh was already promoted");
		}

		// reallocate memory
		Bc = FEM.realloc2dInt(Bc, nb, (p+1)*(p+2)/2);
		T  = FEM.realloc2dInt(T, nt, (p+1)*(p+2)*(p+3)/6);
		// TODO that are too many, exact determination required number of edges
		P  = FEM.realloc2dDouble(P, (nb+nt)*(p+1)*(p+2)/2 + nt*(p-2)*(p-1)*p/6 + nt*2*(p-1), 3 );
		
		// TODO, avoid the use of hashes by using the neigbourhood realtions
		//	(needs edge indices)
		// edge hash 
		Hashtable<Key2,Integer> Eh = new Hashtable<Key2,Integer>(6*nt);
		// face hash
		Hashtable<Key3,Integer> Fh = new Hashtable<Key3,Integer>(4*nt);

		double wp = 1.0/p;
		// permutation arrays
		final int [][] pet = { {0, 1}, {0, 2}, {1, 2} };
		final int [][] pe  = { {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3} };
		final int [][] pf  = { {1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2} };

		// for each boundary face (triangles)
		for (int bdx=0; bdx<nb; bdx++)
		{
			// interior points of the boundary
			int pdx=0;

			// fetch point indices
			int p0 = Bc[bdx][0]-1;
			int p1 = Bc[bdx][1]-1;
			int p2 = Bc[bdx][2]-1;
	
			Key3 key3 = new Key3(p0+1, p1+1, p2+1);
			Fh.put(key3, np);
			//Integer val = Fh.get(key3);

			// face cannot yet have been treated
			// create interior points
			// rows
			for (int i=1; i<p; i++)
			{
				double wr = 1.0/i;

				// for each column
				for (int j=1; j<i; j++)
				{
					P[np][0] = (1-wp*i)*P[p0][0] + (wp*i)*(wr*j*P[p1][0] + (1-wr*j)*P[p2][0]);
					P[np][1] = (1-wp*i)*P[p0][1] + (wp*i)*(wr*j*P[p1][1] + (1-wr*j)*P[p2][1]);
					P[np][2] = (1-wp*i)*P[p0][2] + (wp*i)*(wr*j*P[p1][2] + (1-wr*j)*P[p2][2]);
					np++;
					pdx++;
					Bc[bdx][2 + 3*(p-1) + pdx] = np;
				} // for j columns
			} // for i rows

			// for each edge
			for (int i=0; i<3; i++)
			{
				// prefetch edge indices
				p0 = Bc[bdx][pet[i][0]]-1;
				p1 = Bc[bdx][pet[i][1]]-1;

				// check if edge was already treated
				Key2 key2 = new Key2(p0, p1);
				Integer val = Eh.get(key2);

				if (null == val)
				{
					val = np;
					Eh.put(key2,val);

					// new interior points of this edge
					for (int j=1; j<p; j++)
					{
						P[np][0] = (1-wp*j)*P[p0][0] + (wp*j)*P[p1][0];
						P[np][1] = (1-wp*j)*P[p0][1] + (wp*j)*P[p1][1];
						P[np][2] = (1-wp*j)*P[p0][2] + (wp*j)*P[p1][2];
						np++;
					} // for j (each new point)
				} // if null

				// register point indices
				for (int j=1; j<p; j++)
				{
					Bc[bdx][2 + i*(p-1) + j ] = val + j;
				}
			} // int i (edges)
		} // for each boundary

		// element volumes (tetrahedra)
		for (int tdx=0; tdx<nt; tdx++)
		{
			// edges
			for (int edx=0; edx<6; edx++)
			{
				// prefetch point indices
				int p0 = T[tdx][pe[edx][0]]-1;
				int p1 = T[tdx][pe[edx][1]]-1;

				// edges belong to arbitrarily many tetras, so one can just get, not remove
				Key2 key2 = new Key2(p0, p1);
				Integer val = Eh.get(key2);

				// check of edge was already treatet
				if (null == val)
				{
					val = np;
					Eh.put(key2,val);

					// new interior points of this edge
					for (int j=1; j<p; j++)
					{
						P[np][0] = (1-wp*j)*P[p0][0] + (wp*j)*P[p1][0];
						P[np][1] = (1-wp*j)*P[p0][1] + (wp*j)*P[p1][1];
						P[np][2] = (1-wp*j)*P[p0][2] + (wp*j)*P[p1][2];
						np++;
					} // for j
				} // for j

				// register point indices
				for (int j=1; j<p; j++)
				{
					T[tdx][3 + edx*(p-1) + j ] = val + j;
				} // for j
			} // for edx (edges)

			// faces
			for (int fdx=0; fdx<4; fdx++)
			{
				// prefetch point indices
				int p0 = T[tdx][pf[fdx][0]]-1;
				int p1 = T[tdx][pf[fdx][1]]-1;
				int p2 = T[tdx][pf[fdx][2]]-1;

				// check, wether face has already been treated
				Key3 key3 = new Key3(p0+1, p1+1, p2+1);
				Integer val = Fh.remove(key3);

				if (null == val)
				{
					val = np;
					Fh.put(key3, np);

					// create interior points
					// rows
					for (int i=1; i<p; i++)
					{
						double wr = 1.0/i;
						// columns
						for (int j=1; j<i; j++)
						{
							P[np][0] = (1-wp*i)*P[p0][0] + (wp*i)*(wr*j*P[p1][0] + (1-wr*j)*P[p2][0]);
							P[np][1] = (1-wp*i)*P[p0][1] + (wp*i)*(wr*j*P[p1][1] + (1-wr*j)*P[p2][1]);
							P[np][2] = (1-wp*i)*P[p0][2] + (wp*i)*(wr*j*P[p1][2] + (1-wr*j)*P[p2][2]);
							np++;
						} // for j columns
					} // for i row
				} // if null
				// register the points
				for (int i=1; i<=(p-2)*(p-1)/2; i++)
				{
					T[tdx][3 + 6*(p-1) + fdx*(p-2)*(p-1)/2 + i] = val+i;
				} // for i
			} //for fdx (faces)

			// interior points
			int p0 = T[tdx][0]-1;
			int p1 = T[tdx][1]-1;
			int p2 = T[tdx][2]-1;
			int p3 = T[tdx][3]-1;

			int pdx=0;
			// for each slice
			for (int i=1; i<p; i++)
			{
				double ws = 1.0/i;
				// for each row
				for (int j=1; j<i; j++)
				{
					double wr = 1.0/j;
					// for each column
					for (int k=1; k<j; k++)
					{
						P[np][0] = (1-wp*i)*P[p0][0] + (wp*i)*( (1-ws*j)*P[p1][0] + (ws*j)*((1-wr*k)*P[p2][0] + (wr*k)*P[p3][0]));
						P[np][1] = (1-wp*i)*P[p0][1] + (wp*i)*( (1-ws*j)*P[p1][1] + (ws*j)*((1-wr*k)*P[p2][1] + (wr*k)*P[p3][1]));
						P[np][2] = (1-wp*i)*P[p0][2] + (wp*i)*( (1-ws*j)*P[p1][2] + (ws*j)*((1-wr*k)*P[p2][2] + (wr*k)*P[p3][2]));
						np++;
						pdx++;
						T[tdx][3 + 6*(p-1) + 4*(p-2)*(p-1)/2 + pdx] = np;
					} // for j columns)
				} // for j (rows)
			} // for i (slices)
		} // for tdx (elements)

		// consistency check
		if (!Fh.isEmpty())
		{
			throw new RuntimeException("Mesh_3d.promote: inconsistent mesh");
		}
	} // promote

	public final void demote()
	{
		switch (T[0].length)
		{
			case 10: // second order to first order
				S = new int[8*nt][4];
				ns = 0;
				for (int tdx=0; tdx<nt; tdx++)
				{
					// tetras in corners
					S[ns][0] = T[tdx][0]; S[ns][1] = T[tdx][4]; S[ns][2] = T[tdx][5]; S[ns][3] = T[tdx][6];
					ns++;
					S[ns][0] = T[tdx][1]; S[ns][1] = T[tdx][4]; S[ns][2] = T[tdx][7]; S[ns][3] = T[tdx][8];
					ns++;
					S[ns][0] = T[tdx][2]; S[ns][1] = T[tdx][5]; S[ns][2] = T[tdx][7]; S[ns][3] = T[tdx][9];
					ns++;
					S[ns][0] = T[tdx][3]; S[ns][1] = T[tdx][6]; S[ns][2] = T[tdx][8]; S[ns][3] = T[tdx][9];
					ns++;
					// interior tetras
					S[ns][0] = T[tdx][6]; S[ns][1] = T[tdx][7]; S[ns][2] = T[tdx][8]; S[ns][3] = T[tdx][9];
					ns++;
					S[ns][0] = T[tdx][5]; S[ns][1] = T[tdx][7]; S[ns][2] = T[tdx][6]; S[ns][3] = T[tdx][9];
					ns++;
					S[ns][0] = T[tdx][4]; S[ns][1] = T[tdx][6]; S[ns][2] = T[tdx][7]; S[ns][3] = T[tdx][8];
					ns++;
					S[ns][0] = T[tdx][4]; S[ns][1] = T[tdx][5]; S[ns][2] = T[tdx][7]; S[ns][3] = T[tdx][6];
					ns++;
				} // for tdx
			break;
			case 35: // fourth order to second order
				S = new int[8*nt][10];
				ns = 0;
				for (int tdx=0; tdx<nt; tdx++)
				{
					// tetras in corners
					S[ns][0] = T[tdx][ 0]; S[ns][1] = T[tdx][ 5]; S[ns][2] = T[tdx][ 8]; S[ns][3] = T[tdx][11];
					S[ns][4] = T[tdx][ 4]; S[ns][5] = T[tdx][ 7]; S[ns][6] = T[tdx][10]; S[ns][7] = T[tdx][31];
					S[ns][8] = T[tdx][28]; S[ns][9] = T[tdx][25];
					ns++;
					S[ns][0] = T[tdx][ 1]; S[ns][1] = T[tdx][14]; S[ns][2] = T[tdx][ 5]; S[ns][3] = T[tdx][17];
					S[ns][4] = T[tdx][13]; S[ns][5] = T[tdx][ 6]; S[ns][6] = T[tdx][16]; S[ns][7] = T[tdx][33];
					S[ns][8] = T[tdx][22]; S[ns][9] = T[tdx][30];
					ns++;
					S[ns][0] = T[tdx][ 2]; S[ns][1] = T[tdx][ 8]; S[ns][2] = T[tdx][14]; S[ns][3] = T[tdx][20];
					S[ns][4] = T[tdx][ 9]; S[ns][5] = T[tdx][15]; S[ns][6] = T[tdx][19]; S[ns][7] = T[tdx][32];
					S[ns][8] = T[tdx][27]; S[ns][9] = T[tdx][24];
					ns++;
					S[ns][0] = T[tdx][ 3]; S[ns][1] = T[tdx][11]; S[ns][2] = T[tdx][17]; S[ns][3] = T[tdx][20];
					S[ns][4] = T[tdx][12]; S[ns][5] = T[tdx][18]; S[ns][6] = T[tdx][21]; S[ns][7] = T[tdx][29];
					S[ns][8] = T[tdx][26]; S[ns][9] = T[tdx][23];
					ns++;
					// interior tetra
					S[ns][0] = T[tdx][ 5]; S[ns][1] = T[tdx][ 8]; S[ns][2] = T[tdx][14]; S[ns][3] = T[tdx][20];
					S[ns][4] = T[tdx][31]; S[ns][5] = T[tdx][33]; S[ns][6] = T[tdx][27]; S[ns][7] = T[tdx][32];
					S[ns][8] = T[tdx][34]; S[ns][9] = T[tdx][24];
					ns++;
					S[ns][0] = T[tdx][ 8]; S[ns][1] = T[tdx][11]; S[ns][2] = T[tdx][20]; S[ns][3] = T[tdx][ 5];
					S[ns][4] = T[tdx][25]; S[ns][5] = T[tdx][27]; S[ns][6] = T[tdx][31]; S[ns][7] = T[tdx][26];
					S[ns][8] = T[tdx][28]; S[ns][9] = T[tdx][34];
					ns++;
					S[ns][0] = T[tdx][ 5]; S[ns][1] = T[tdx][17]; S[ns][2] = T[tdx][11]; S[ns][3] = T[tdx][20];
					S[ns][4] = T[tdx][30]; S[ns][5] = T[tdx][28]; S[ns][6] = T[tdx][34]; S[ns][7] = T[tdx][29];
					S[ns][8] = T[tdx][23]; S[ns][9] = T[tdx][26];
					ns++;
					S[ns][0] = T[tdx][14]; S[ns][1] = T[tdx][20]; S[ns][2] = T[tdx][17]; S[ns][3] = T[tdx][ 5];
					S[ns][4] = T[tdx][24]; S[ns][5] = T[tdx][22]; S[ns][6] = T[tdx][33]; S[ns][7] = T[tdx][23];
					S[ns][8] = T[tdx][34]; S[ns][9] = T[tdx][30];
					ns++;
				} // for tdx
			break;
			default:
				throw new RuntimeException("Not yet implemented or impossible");
		} // switch T[0].length
	} // demote
} // class Mesh_3D

