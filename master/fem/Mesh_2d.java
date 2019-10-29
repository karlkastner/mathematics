// Sat Jun 16 16:31:54 MSK 2012
// Karl Kästner, Berlin
import Jama.*;
import java.util.Hashtable;

// TODO doxygen comments
// TODO counters np, nt, nb

public final class Mesh_2d
{
	// point coordinates
	public double [][] P;
	public int np;
	// elements (triangles)
	public int [][] T;
	public int nt;
	// alternative triagulation with basis function p/2 for error estimation
	public int [][] S;
	public int ns;
	// boundary sides
	public int [][] Bc;
	public int nb;
	// neighbour indices
	public int [][] N;
	// dimension
	final int DIM = 2;

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

	// default constructor
	public Mesh_2d() {};

	public Mesh_2d(double [][] P, int [][] T, int [][] Bc)
	{
		this.P = P;
		this.np = P.length;
		this.T = T;
		this.nt = T.length;
		this.Bc = Bc;
		this.nb = Bc.length;
		// TODO, null other arrays
		this.P_local = null;
		this.Phi = null;
	} // constructor

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

		// there are exactly p(p-1)/2 new points for each triangle
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
			Matrix mC = new Matrix(Phi[idx]);
			
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
} // class Mesh

