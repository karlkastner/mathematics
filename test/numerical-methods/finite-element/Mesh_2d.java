// Sat Jun 16 16:31:54 MSK 2012
// Karl Kästner, Berlin
import Jama.*;
import java.util.Hashtable;
import java.util.*;
import java.lang.*;

// TODO doxygen comments
// TODO counters np, nt, nb
// TODO naming convention, this is just for triangles

public final class Mesh_2d extends Mesh
{
	// point coordinates
	// TODO either split this up to X Y Z or make first dimension the smaller one
//	public double [][] P;
//	public int np;
	// elements (triangles)
//	public int [][] T;
//	public int nt;
	// alternative triagulation with basis function p/2 for error estimation
	public int [][] S;
	public int ns;
	// boundary sides
//	public int [][] Bc;
//	public int nb;
	// neighbour indices
//	public int [][] N;
	// dimension
//	final int DIM = 2;
	// properties derived from primery properties
	// test function coefficients
	public double [][][] Phi;
	// prefetched points per element
	public double [][][] P_local;
	// side length per element
	public double h_side[];
	// interior angles per element
	public double degen[];
	// determinant of each element
	public double determinant[];
	// centre coordinates of each element
	public double C[][];

	public double l_boundary[];

	// tree of bubbles to assign random points to the containing triangle
	// in logarithmic time (for interpolation)
	Qtree qtree;

	// default constructor
	public Mesh_2d()
	{
		DIM = 2;
	};

	public Mesh_2d(final int np, final int nt, final int nbc)
	{
		P  = new double [np][2];
		T  = new int [nt][3];
		Bc = new int [nbc][2];
	}

	public Mesh_2d(double [][] P, int [][] T, int [][] Bc)
	{
		DIM = 2;
		this.P = P;
		if (null == P)
		{
			np = 0;
		} else {
			np = P.length;
		}
		this.T = T;
		if (null == T)
		{
			nt = 0;
		} else {
			nt = T.length;
		}
		this.Bc = Bc;
		if (null == Bc)
		{
			this.nb = 0;
		} else {
			this.nb = Bc.length;
		}
		// TODO, null other arrays
		this.P_local = null;
		this.Phi     = null;
	} // constructor

	Mesh_2d(final double [][] P)
	{
		this.P = P; // legal ?
		// initial delaunay
		delaunay();
		// improvement
	}

	public final int add_T(final int p0, final int p1, final int p2)
	{
		if (p0 == p1 || p0 == p2 || p1 == p2) throw new RuntimeException("Error identical point indices");
		T = AArray.resizeI(T,nt+1);
		T[nt][0] = p0;
		T[nt][1] = p1;
		T[nt][2] = p2;
		nt++;
		return nt;
	} // add_T()

	public final void add_Bc(final int p0, final int p1)
	{
		if (p0 == p1) throw new RuntimeException();
		Bc = AArray.resizeI(Bc,nb+1);
		Bc[nb][0] = p0;
		Bc[nb][1] = p1;
		nb++;
	} // add_Bc()

	// TODO make N and element_neighbours part of the Mesh_2d class
	// compute neighbouring elements
	public final void element_neighbours()
	{
		int lt1 = nt; //T.length;
		int lb1 = nb; //Bc.length;

		// allocate memory
		N = new int[lt1][DIM+1];

		// hash of triangle sides
		// key: lower_point_index, higher_point_index
		// value: triangle index, index of opposit point
		// if boundary, boundary index and negative number of boundary
		Hashtable<Key2,int[]> Sh = new Hashtable<Key2,int[]>((DIM+1)*lt1);

		// push the boundaries
		for (int idx=0; idx<lb1; idx++)
		{
			// push row and column in N, negative, domain
			int val[] = new int[1];
			val[0] = -idx-2;
			Sh.put(new Key2(Bc[idx][0], Bc[idx][1]), val);
		} // for idx

		for (int idx=0; idx<lt1; idx++)
		{
			for (int jdx=0; jdx<3; jdx++)
			{
				// point at opposit side
				Key2 key = new Key2(T[idx][(jdx+1) % 3], T[idx][(jdx+2) % 3]);
				int [] val = Sh.remove(key);
				if (null == val)
				{
					val = new int[2];
					val[0] = idx; val[1] = jdx;
					Sh.put(key, val);
				} else {
					// Neighbourhood relation found
					N[idx][jdx] = val[0]+1;
					if (val[0] >= 0)
					{
						N[val[0]][val[1]] = idx+1;
					}
				}
			} // for jdx
		} // for idx
		
		if (!Sh.isEmpty())
		{
			// error('fem_2d_element_boundary','inconsistent mesh');
			System.err.println("element_neighbours: inconsistent mesh");
			//throw new Throwable();
		}
	} // element_neighbours


	public final void prefetch()
	{
		int nt1 = nt; //T.length;
		int nt2 = T[0].length; // dangerous if T empty
		// determine number of the power (result of expression is alway integer)
		int nv = (int) (-1.5 + Math.sqrt(2.0*nt2 + 0.25));

		// allocate memory
		Phi         = new double[nt1][nt2][nt2];
//		h_side      = new double[nt1][DIM+1];
//		s_angle     = new double[nt1][DIM+1];
		h_side      = new double[nt1];
		degen       = new double[nt1];
		C           = new double[nt1][DIM];
		determinant = new double[nt1];
//		l_boundary  = new double[nb]; //Bc.length];
		P_local     = new double[nt1][DIM+1][DIM];
		
		// temporary arrays
		double [][]   A    = new double[DIM+1][DIM+1];
		double [][]   A_   = new double[nt2][DIM];
		double [][]   Va   = new double[nt2][(nv+1)*(nv+2)/2];
		Matrix        mVa  = new Matrix(Va);
		Matrix        mA   = new Matrix(A);
		
		for (int tdx=0; tdx<nt1; tdx++)
		{
			// TODO, use same shortcut in other functions
			int [] T_tdx = T[tdx];

			// prefetched point coordinates (TODO, use below and do not rely on compiler)
			P_local[tdx][0][0] = P[T_tdx[0]-1][0]; P_local[tdx][0][1] = P[T_tdx[0]-1][1];
			P_local[tdx][1][0] = P[T_tdx[1]-1][0]; P_local[tdx][1][1] = P[T_tdx[1]-1][1];
			P_local[tdx][2][0] = P[T_tdx[2]-1][0]; P_local[tdx][2][1] = P[T_tdx[2]-1][1];

			//  prefetch triangle points
			for (int j=0; j<nt2; j++)
			{
				A_[j][0] = P[T_tdx[j]-1][0];
				A_[j][1] = P[T_tdx[j]-1][1];
			}
			
			// construct the Vandermonde matrix
			FEM.vander_2d(Va, A_, nv);
			
			// invert the Vandermonde matrix
			// TODO fetch singular matrices and use SVD to compute the pseudoinverse
			Matrix mC = mVa.inverse();

			// write to output array
			Phi[tdx] = mC.getArray();

			// calculate element centre coordinate
			C[tdx][0] = 1.0/3.0*(P[T_tdx[0]-1][0] + P[T_tdx[1]-1][0] + P[T_tdx[2]-1][0]);
			C[tdx][1] = 1.0/3.0*(P[T_tdx[0]-1][1] + P[T_tdx[1]-1][1] + P[T_tdx[2]-1][1]);

			// fetch corner points
			// TODO do not fecht again, use submatrix
			A[0][0] = 1; A[0][1] = P[T_tdx[0]-1][0]; A[0][2] = P[T_tdx[0]-1][1];
	                A[1][0] = 1; A[1][1] = P[T_tdx[1]-1][0]; A[1][2] = P[T_tdx[1]-1][1];
	                A[2][0] = 1; A[2][1] = P[T_tdx[2]-1][0]; A[2][2] = P[T_tdx[2]-1][1];

			// calculate element determinant
			determinant[tdx] = Math.abs(mA.det());

			// todo: warn for clockwise (area<0) and degenerated triangles (area/h_max < eps)

			// length of side opposit the points
/*
			double a = Math.sqrt( (P[T[tdx][1]-1][0] - P[T[tdx][2]-1][0])*(P[T[tdx][1]-1][0] - P[T[tdx][2]-1][0]) 
						+ (P[T[tdx][1]-1][1] - P[T[tdx][2]-1][1])*(P[T[tdx][1]-1][1] - P[T[tdx][2]-1][1]) );
			double b = Math.sqrt((P[T[tdx][0]-1][0] - P[T[tdx][2]-1][0])*(P[T[tdx][0]-1][0] - P[T[tdx][2]-1][0])
						+ (P[T[tdx][0]-1][1] - P[T[tdx][2]-1][1])*(P[T[tdx][0]-1][1] - P[T[tdx][2]-1][1]));
			double c = Math.sqrt((P[T[tdx][0]-1][0] - P[T[tdx][1]-1][0])*(P[T[tdx][0]-1][0] - P[T[tdx][1]-1][0])
						+ (P[T[tdx][0]-1][1] - P[T[tdx][1]-1][1])*(P[T[tdx][0]-1][1] - P[T[tdx][1]-1][1]));
*/
			//TODO as in 3d, a one colum measure for h and degeneracy should be computed her
			// Geomean is probably not a good idea
			double h = 0;
			double a = FEM.dist2(P[T_tdx[0]-1], P[T_tdx[1]-1]);
			double b = FEM.dist2(P[T_tdx[0]-1], P[T_tdx[2]-1]);
			double c = FEM.dist2(P[T_tdx[1]-1], P[T_tdx[2]-1]);
			double [] c_ = new double[2];
			// todo, no magick numbers
			switch (1)
			{
				case 0 :
					// max norm
					h =             a;  //FEM.dist2(P[T_tdx[0]-1], P[T_tdx[1]-1]);
					h = Math.max(h, b); //FEM.dist2(P[T_tdx[0]-1], P[T_tdx[1]-1]));
					h = Math.max(h, c); //FEM.dist2(P[T_tdx[0]-1], P[T_tdx[1]-1]));
				break;
				case 1: 
					// diameter of the circumcircle
					c_[0] = 1.0/3.0*( P[T_tdx[0]-1][0] + P[T_tdx[1]-1][0] + P[T_tdx[2]-1][0]); 
					c_[1] = 1.0/3.0*( P[T_tdx[0]-1][1] + P[T_tdx[1]-1][1] + P[T_tdx[2]-1][1]);
					h = 2*FEM.dist2(c_, P[T_tdx[0]-1]);
				break;
			}
			h_side[tdx] = h;

			// degeneracy
			double sigma = 1;
			switch (1)
			{
				case 0:
					// sin_alpha = 0.5|A|/(0.5*a*b)
					sigma =                 a*b;
					sigma = Math.max(sigma, b*c);
					sigma = Math.max(sigma, a*c);
				break;
				case 1:	// area / diameter_circumsphere^2
					sigma = h*h;
				break;
			}
			degen[tdx] = determinant[tdx]/sigma;

			/*

			h_side[tdx][0] = a;
			h_side[tdx][1] = b;
			h_side[tdx][2] = c;

			// sine of interior angles at the points
			double sin_a = 0.5*determinant[tdx]/(0.5*b*c);
			double sin_b = 0.5*determinant[tdx]/(0.5*a*c);
			double sin_c = 0.5*determinant[tdx]/(0.5*a*b);

			s_angle[tdx][0] = sin_a;
			s_angle[tdx][1] = sin_b;
			s_angle[tdx][2] = sin_c;
			*/
		} // for i

		// boundary length
		// TODO this is actually not necessary to compute, make it a test case
		/*
		for (int idx=0; idx<Bc.length; idx++)
		{
			l_boundary[idx] = Math.sqrt( (P[Bc[idx][0]-1][0] - P[Bc[idx][1]-1][0])*(P[Bc[idx][0]-1][0] - P[Bc[idx][1]-1][0])
						   + (P[Bc[idx][0]-1][1] - P[Bc[idx][1]-1][1])*(P[Bc[idx][0]-1][1] - P[Bc[idx][1]-1][1]) );
		}
		*/
	} // prefetch_2d
	
	// Promote
	// Sat Jun 16 16:56:14 MSK 2012
	// Karl Kästner, Berlin

	// add points required for Lagrangian basis functions of order p
	public final void promote(final int p)
	{
/*
		// get array dimensions
		int np = P.length;
		int nt = T.length;
		int nb = Bc.length;
*/

		// preallocate memory
		// each boundary element has p+1 points
		Bc = FEM.realloc2dInt(Bc, nb, p+1);

		// there is exactly p(p-1)/2 new points for each triangle
		// and p-1 new points for each boundary segment minus 1
		P = FEM.realloc2dDouble(P, ((p*p-1)*nt + (p-1)*nb)/2 + np, DIM); //P[0].length);
//		P = FEM.realloc2dDouble(P, np + nt*p*(p+1)/2 + nb*(p-1) - 1, P[0].length);

		// each triangle hash (p+1)*(p+2)/2 points
		T = FEM.realloc2dInt(T, nt, (p+1)*(p+2)/2);
	
		// triangle egde hashtable
		// there will never be more than np/2 entries in the hashtable
		// TODO, use neighbourhood relation and not a hash
		Hashtable<Key2,Integer> H = new Hashtable<Key2,Integer>(np);

		double wp = 1.0/p;

		// push the boundaries points
		for (int idx=0; idx<nb; idx++)
		{
			// fetch point indices
			int p0 = Bc[idx][0]-1;
			int p1 = Bc[idx][1]-1;

			// for each new point on the boundary segment
			for ( int i=1; i<p; i++)
			{
				P[np+i-1][0] = i*wp*P[p0][0] + (p-i)*wp*P[p1][0];
				P[np+i-1][1] = i*wp*P[p0][1] + (p-i)*wp*P[p1][1];
			} // for each new point

			// add new point indices
			for (int i=1; i<p; i++)
			{
				Bc[idx][1 + i] = np + i;
			} // for each new point

			// memorise the new segment
			H.put(new Key2(p0,p1),np);
			np += p-1;
		} // for idx

		int [][] f = { { 0, 1 }, {0, 2}, {1, 2} };

		// for each triangle
		for (int tdx=0; tdx<nt; tdx++)
		{
			// interior points
			int pdx = 0;

			// fetch point indices
			int p0 = T[tdx][0]-1;
			int p1 = T[tdx][1]-1;
			int p2 = T[tdx][2]-1;

			// for each row
			for (int i=1; i<p; i++)
			{
				double wc = 1.0/i;

				// for each column
				for (int j=1; j<i; j++)
				{
					// new point coordinate
					P[np][0] = (1-wp*i)*P[p0][0] + wp*i*(wc*j*P[p1][0] + (1-wc*j)*P[p2][0]);
					P[np][1] = (1-wp*i)*P[p0][1] + wp*i*(wc*j*P[p1][1] + (1-wc*j)*P[p2][1]);
					np++;
					pdx++;
					T[tdx][2 + 3*(p-1) + pdx] = np;
				} // for j
			} // for j

			// for each edge
			for (int edx=0; edx<3; edx++)
			{
				// get end points of the edge
				p0 = T[tdx][f[edx][0]]-1;
				p1 = T[tdx][f[edx][1]]-1;

				// check wether the points exist already
				Key2 key2 = new Key2(p0, p1);
				Integer val = H.remove(key2);
			
				if (null == val)
				{
					// create the interior points, as they do not yet exist
					val = np;
					H.put(key2, val);
				
					// for each new point on the edge
					for ( int i=1; i<p; i++)
					{
						P[np][0] = i*wp*P[p0][0] + (p-i)*wp*P[p1][0];
						P[np][1] = i*wp*P[p0][1] + (p-i)*wp*P[p1][1];
						np++;
					} // for each new point
				} // if null
			
				// add new point indices
				for (int i=1; i<p; i++)
				{
					T[tdx][2 + edx*(p-1) + i ] = val + i;
				} // for each new point
			} // for edx, each edge
		} // for tdx, each triangle

		// verify edge consistency
		if (!H.isEmpty())
		{
			throw new RuntimeException("promote_3_21: inconsitent mesh");
		}
	} // promote

	public final void demote()
	{
//		int nt = T.length;
		switch (T[0].length)
		{
			case 6:	// second order to first order
				S = new int[4*nt][3];
				ns = 0;
				for (int tdx=0; tdx<nt;tdx++)
				{
					// triangles in corners
					S[ns][0] = T[tdx][0]; S[ns][1] = T[tdx][3]; S[ns][2] = T[tdx][4];
					ns++;
					S[ns][0] = T[tdx][1]; S[ns][1] = T[tdx][3]; S[ns][2] = T[tdx][5];
					ns++;
					S[ns][0] = T[tdx][2]; S[ns][1] = T[tdx][4]; S[ns][2] = T[tdx][5];
					ns++;
					// interior triangle
					S[ns][0] = T[tdx][3]; S[ns][1] = T[tdx][4]; S[ns][2] = T[tdx][5];
					ns++;
				} // for tdx 
			break;
			case 15: // fourth order to second order
				S = new int[4*nt][6];
				ns = 0;
				for (int tdx=0; tdx<nt;tdx++)
				{
					// triangles in corners
					S[ns][0] = T[tdx][0]; S[ns][1] = T[tdx][4]; S[ns][2] = T[tdx][7];
					S[ns][3] = T[tdx][3]; S[ns][4] = T[tdx][6]; S[ns][5] = T[tdx][12];
					ns++;
					S[ns][0] = T[tdx][1]; S[ns][1] = T[tdx][10]; S[ns][2] = T[tdx][4];
					S[ns][3] = T[tdx][9]; S[ns][4] = T[tdx][ 5]; S[ns][5] = T[tdx][13];
					ns++;
					S[ns][0] = T[tdx][2]; S[ns][1] = T[tdx][ 7]; S[ns][2] = T[tdx][10];
					S[ns][3] = T[tdx][8]; S[ns][4] = T[tdx][11]; S[ns][5] = T[tdx][14];
					ns++;
					// interior triangle
					S[ns][0] = T[tdx][ 4]; S[ns][1] = T[tdx][10]; S[ns][2] = T[tdx][ 7];
					S[ns][3] = T[tdx][13]; S[ns][4] = T[tdx][12]; S[ns][5] = T[tdx][14];
					ns++;
				} // for tdx
			break;
/*
			case 36: // sith order to third order
			break;
*/
			default:
				throw new RuntimeException("Not yet implemented or impossible");
		} // switch T[0].length
	} // demote

	// computes the highest (d-th) derivative
	// TODO, check for prefetch here
	public final double [][] dV(final double [] V, final int d)
	{
		int lt1 = nt; //T.length;

		// allocate memory
		int n = (d+1)*(d+2)/2;
		double [][] dV      = new double[lt1][d+1];
		double [][] V_local = new double[n][1];
		Matrix mV_local     = new Matrix(V_local);

		// compute the coefficients
		final int p = d;
        	int [] f = new int[p+1];
		f[0] = 1;
        	for (int pdx=1; pdx<=p; pdx++)
		{
                	f[pdx] = pdx*f[pdx-1];
        	}

        	final int n0 = p*(p+1)/2;
		int [] dC = new int[p+1];

		int ndx=0;
                for (int dx=p; dx>=0; dx--)
		{
                       	int dy=p-dx;
                        dC[ndx] = f[dx]*f[dy];
                        ndx++;
		} // dx

		// for each element
		for (int idx=0; idx<lt1; idx++)
		{
			// fetch the function values
			for (int jdx=0; jdx<n; jdx++)
			{
				V_local[jdx][0] = V[T[idx][jdx]-1];
			}
			// test/trial function coefficients
			Matrix mC;
			if (null != Phi)
			{
				mC = new Matrix(Phi[idx]);
			} else {
				
				int nt2 = T[0].length; //TODO dangerous if T is empty
				int nv = (int) Math.round(-1.5 + Math.sqrt(2*nt2 + 0.25));
				double[][] A_   = new double[nt2][2];
				double[][] Va   = new double[nt2][(nv+1)*(nv+2)/2];
				Matrix mVa   = new Matrix(Va);
				// fetch all point coordinates, not just corner point coordinates
				// TODO, these coordinates can be computed on the fly (except on curved boundaries)
				for (int j=0; j<nt2; j++)
				{
					A_[j][0] = P[T[idx][j]-1][0];
					A_[j][1] = P[T[idx][j]-1][1];
				}
				// construct the triangle point Vandermonde matrix spanning the test function polynomials
				// structure: 1 x y xy x^2 y^2 x^2y xy^2 x^3 y^3
				FEM.vander_2d(Va, A_, nv);
				//FEM.vander_2d(Va, mA.getMatrix(0, lt2-1, 1, 2).getArray(), nv);
	
				// compute test functions
				// C(i,:) * A(:,i) = 1; C(i,:) * A(:,j) = 0, i<>j
				mC = mVa.inverse();
				//C = mC.getArray();
			}
			
			Matrix mD = mC.times(mV_local);
			double [][] D = mD.getArray();

			// compute the coefficients of the partial derivatives
			for (int jdx=0; jdx<dC.length; jdx++)
			{
				dV[idx][jdx] = dC[jdx]*D[n0+jdx][0];
			}

			/*
			// compute dth order partial derivative
			switch (d)
			{
				case 1:
					// extract coefficients of first order partial derivatives
					dV[idx][0] = D[1][0];
					dV[idx][1] = D[2][0];
				break;
				case 2:
					// extract coefficients of second order partial derivatives
					dV[idx][0] = 2*D[3][0];
					dV[idx][1] =   D[4][0];
					dV[idx][2] = 2*D[5][0];
				break;
				case 3:
					// extract coefficients of third order partial derivatives
					dV[idx][0] = 6*D[6][0];
					dV[idx][1] = 2*D[7][0];
					dV[idx][2] = 2*D[8][0];
					dV[idx][3] = 6*D[9][0];
				break;
				case 4:
					// extract coefficients of fourth order partial derivatives
					dV[idx][0] = 24*D[10][0];
					dV[idx][1] =  6*D[11][0];
					dV[idx][2] =  4*D[12][0];
					dV[idx][3] =  6*D[13][0];
					dV[idx][4] = 24*D[14][0];
				break;
				case 5:
					// extract coefficients of fith order partial derivatives
					dV[idx][0] = 120*D[15][0];
					dV[idx][1] =  24*D[16][0];
					dV[idx][2] =  12*D[17][0];
					dV[idx][3] =  12*D[18][0];
					dV[idx][4] =  24*D[19][0];
					dV[idx][5] = 120*D[20][0];
				break;
				default:
					throw new RuntimeException("Not yet implemented");
			} // switch (d)
			*/
		} // for int i

		return dV;
	} // dV

	// estimates the error
	public final Object estimate_error(final  double[][] dV, final int m, final int k, final int s)
	{
		//  get length
		int lt1 = nt; //T.length;

		//  allocate memory
		double err_max = 0;
		double thresh  = 0;
		double err []  = new double[lt1];
		double nH  []  = new double[lt1];
		Object obj []  = new Object[4];

		// observe convergence of the sth-order derivative
		// see Strang, Fix Theorem 3.7, An Analysis of the Finite Element Method
		double rate = Math.min(k-s, 2*(k-m));

		// constants
		// C_ = [1/3 1/38 1/870 1/35430 1/882847]; % for poisson mode 0
		// C_ = [1/3, 1/13, 1/87, 1/908, 1/12360]; % for Schrödinginger mode 0
		//double [] C_ = { 1.0/9, 1.0/(3*13), 1.0/(3*87), 1.0/(3*908), 1.0/(3*12360) }; // % for Schrödinger mode 1
		double [] C_ = { 1.0/3.0, 1.0/13.0, 1.0/87.0, 1.0/908.0, 1.0/12360.0 }; // % for Schrödinger mode 1

		// prefetch
		if (null == C)
		{
			C = new double[lt1][2];
			for (int tdx=0; tdx<lt1; tdx++)
			{
				int [] T_tdx = T[tdx];
				// calculate element centre coordinate
				C[tdx][0] = 1.0/3.0*(P[T_tdx[0]-1][0] + P[T_tdx[1]-1][0] + P[T_tdx[2]-1][0]);
				C[tdx][1] = 1.0/3.0*(P[T_tdx[0]-1][1] + P[T_tdx[1]-1][1] + P[T_tdx[2]-1][1]);
			}
		}

		// calculate estimated norm of the second derivative per element
		for (int idx=0; idx<lt1; idx++)
		{
			// seminorm of derivativ p+1
			// for all three triangle neighbours
			// Extension to Eriksson, Johnson 1988, Adaptive Finite Element Method for Linear Elliptic Problems
			// fails for s=0 with nonsmooth data
			for (int jdx=0; jdx<3; jdx++)
			{
				// check that neighbour exists (not a domain boundary)
				if (N[idx][jdx] > 0)
				{
					double dx      = C[idx][0] - C[N[idx][jdx]-1][0];
					double dy      = C[idx][1] - C[N[idx][jdx]-1][1];
					double dr_sqr  = dx*dx + dy*dy;
					double nH_ = 0;
					// take the maximum of the partial derivatives
					for (int ddx=0; ddx<dV[0].length; ddx++)
					{
						nH_ = Math.max(nH_,Math.abs((dV[idx][ddx] - dV[N[idx][jdx]-1][ddx])*dx)/dr_sqr);
						nH_ = Math.max(nH_,Math.abs((dV[idx][ddx] - dV[N[idx][jdx]-1][ddx])*dy)/dr_sqr);
					}
					// take maximum maximum over of all three neighbours
					nH[idx] = Math.max(nH[idx], nH_);
				} // N(idx,jdx) > 0
			} // for jdx

			// local error of the element
			err[idx] = C_[k-2] * nH[idx] * Math.pow( h_side[idx] / degen[idx], rate);
//				* Math.pow(h_side[idx][0]*h_side[idx][1]*h_side[idx][2]
//						/ (s_angle[idx][0]*s_angle[idx][1]*s_angle[idx][2]), 1.0/3.0*rate);

			// global error over all elements
			err_max = Math.max(err_max, err[idx]);
		} // for idx

		// threshold for refinement
		thresh = Math.pow(0.5,rate)*err_max;

		obj[0] = err;
		obj[1] = err_max;
		obj[2] = thresh;
		obj[3] = nH;	

		return obj;
	} // estimate_error

	// delaunay triangulation in O(log(n)n) time
	private void delaunay()
	{
		int flipmode = 0;
		int insmode = 0;

		// TODO neighbourhood relations

		// check that there are at least three points
		if (P.length < 3)
		{
			throw new RuntimeException("Set of input points must contain at least of size three.");
		}
	
		// form the first triangle out of the first three points
		add_T(1,2,3);
	
		// prepare point indices
		int Q[] = new int[P.length-3];
		for (int idx=0; idx<P.length; idx++)
		{
			Q[idx] = idx+3;
		} // for idx
	
		// assign all points to the first triangle
		// TODO, check, that all points are contained
		// one could do that by forming a triangle around the circumcircle first
		// push first parent triangle containing all points
		java.util.Stack<int[]> T_stack = new java.util.Stack<int[]>();
		if (P.length-3 > 0)
		{
			int [] T_top = new int[3];
			T_top[0] = 0; // index of first triangle in triangle array
			T_top[1] = 3; // index of first contained point in point array
			T_top[2] = P.length; // index of last+1 contained point in point array
			T_stack.push(T_top);
		}
	
		if (1 == insmode)
		{
			// TODO find point closest to centre of mother triangle
		}
	
		// while there are unsplit triangles containing unconnected points
		while (!T_stack.empty())
		{
			// get next triangle index
			int [] T_top = T_stack.pop();
			int tdx = T_top[0];
			int p_top = T_top[1];
	
			// split triangle into three subtriangles
			// TODO, add children pointer or remove parent triangle
			int c0 = add_T(T[tdx][0], T[tdx][1],  Q[p_top]);
			int c1 = add_T(T[tdx][0],   Q[p_top], T[tdx][2]);
			int c2 = add_T( Q[p_top], T[tdx][1], T[tdx][2]);
	
			// find points contained in the first child triangle
			int np = 0;
			for (int pdx=p_top; pdx<T_top[2]; pdx++)
			{
				if (contains(c0, P[pdx])); //Contains.contains2D(P,T,c0,pdx))
				{
					AArray.swapI(Q,pdx,p_top+np);
					np++;
				}
			}
	
			// push the first triangle onto the stack if it contains2D any points 
			if (np != 0)
			{
				int [] T_next = new int[3];
				T_next[0] = c0;
				T_next[1] = p_top;
				T_next[2] = p_top+np;
				T_stack.push(T_next);
			}
			
			// find points contained in second child triagnle
			p_top += np;
			np = 0;
			for (int pdx=p_top; pdx<T_top[2]; pdx++)
			{
				if (contains(c1,P[pdx])) // Contains.contains2D(P,T,c1,pdx))
				{
					AArray.swapI(Q, pdx, p_top+np);
					np++;
				}
			}
	
			if (np != 0)
			{
				int [] T_next = new int[3];
				T_next[0] = c1;
				T_next[1] = p_top;
				T_next[2] = p_top+np;
				T_stack.push(T_next);
			}
	
			// remaining points are contained in the third triangle
			p_top += np;
			if (p_top != T_top[2])
			{
				int [] T_next = new int[3];
				T_next[0] = c2;
				T_next[1] = p_top;
				T_next[2] = T_top[2];
				T_stack.push(T_next);
			}
		
			// flip sides of this triangle	
			if (0 == flipmode)
			{
				flip_recursively(P,T,N,c0);
				flip_recursively(P,T,N,c1);
				flip_recursively(P,T,N,c2);
			}
	
		} // triangle splitting recursion
				
		// flip triangles afterwards
		if (1 == flipmode)
		{
			// for each triangle
			for (int tdx = 0; tdx < nt; tdx++)
			{
				// flip recursively
				flip_recursively(P, T, N, tdx);
			}	
		} // if 1 == flipmode
	} // void delaunay

	void flip_recursively(final double [][] P, int [][] T, int [][] N, final int tdx)
	{
		// try to flip edge 1
		if (false)
		{
			// flip the connecting side
			// update neighbourhood relations
			// recursively flip both triangles
		/*	flip_recursively(P,T,N,n0);
			flip_recursively(P,T,N,n1);
			flip_recursively(P,T,N,n2);
			flip_recursively(P,T,N,m0);
			flip_recursively(P,T,N,m1);
			flip_recursively(P,T,N,m2);
		*/
		} else if (false) { // why else ?
			// try to flip edge 2
			// TODO
		} else if (false) {
		// try to flip with neighbour 3
			// TODO
		}
	} // void flip()
	
	// Sat Mar 15 14:55:01 WIB 2014
	// TODO this could be optimised by splitting the domain recursively in form of a bounding box tree
	// current run time nt*np
	// tree run time log(n)*n
	public ArrayList [] assign_points(double [] X, double [] Y)
	{
		ArrayList [] assoc = new ArrayList[nt];
		for (int idx=0; idx<nt; idx++)
		{
			assoc[idx] = new ArrayList<Integer>();
		}
		qtree = new Qtree(this);
		//qtree.plot();

		// for each point
		int m = 0;
		for (int idx=0; idx<X.length; idx++)
		{
			double [] p = { X[idx], Y[idx] };
			int tdx = qtree.contains(p);
			if (tdx > 0)
			{
				assoc[tdx].add(new Integer(idx+1));
				m++;
			}
		}
		return assoc;
	} // assign_points

	public int contains(final double [] p)
	{
		if (null == qtree)
		{
			qtree = new Qtree(this);
		}
		int retval = qtree.contains(p)+1;
		return retval;
	}

	public boolean contains(final int tdx, final double [] p)
	{
		double A[][] = { { 1.0, 1.0, 1.0 },
       	                         { P[T[tdx][0]-1][0], P[T[tdx][1]-1][0], P[T[tdx][2]-1][0] },
       	                         { P[T[tdx][0]-1][1], P[T[tdx][1]-1][1], P[T[tdx][2]-1][1] } };
		double [][] b = {{1.0}, {p[0]}, {p[1]}};
		Matrix mA      = new Matrix(A);
		Matrix mP      = new Matrix(b);
		Matrix mC      = mA.solve(mP);
		double [][] c  = mC.getArray();	

		boolean retval = (c[0][0] >= 0.0) && (c[1][0] >= 0.0) && (c[2][0] >= 0.0)
		                 && (c[0][0] <= 1.0) && (c[1][0] <= 1.0) && (c[2][0] <= 1.0);
//		System.out.println(c[0][0] + " " + c[1][0] + " " + c[2][0] + " " + retval );

		return retval;	
	} // contains2D
} // class Mesh

