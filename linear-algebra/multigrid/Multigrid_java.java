// Tue 23 Jan 20:40:10 CET 2024
// javac -target 1.8 Multigrid.java 
import java.util.Arrays;  
import Jama.Matrix;
//
// solve reaction-advection-diffusion system on a square nxn grid in two dimensions
// n must be a power of 2
//
//	A*z + (adx*Dx + ady*Dy + dx*Dx^2 + dy*Dy^2)*z = b
//
// values of z are stored in a flat 1xn^2 vector
//
// x can consist of multiple (nvar) coupled variables with each nxn values:
//
//	[A11 ... A1n][ x1]   [d1*D^2        ][ x1]   [ b1]
//	[    ...    ][...] + [    ...       ][...] = [...]
//	[An1 ... Ann][ xn]   [        dn*D^2][ xn]   [ bn]
//
// reaction coefficients Aij can be constant 1x1 or variable (nxn)
// diffusion coefficients di have to be constant but can differ between the variables xi
//
// the reaction coefficients Aij are resampled algebraically, 
// the diffusion part D2 is resampled geometrically
//
// Example solving a single backward-euler time step of the heat equation:
//	 dx/dt = A_*x - d*D2*x
//	 (I - dt*A_ + dt*d_*D2) x_kt+1 = x_kt
// Therefore Aii = 1 - dt*A_ii, and d = dt*dii
//
// TODO speed up by processing 2D matvecs in blocs
class Multigrid_java {
		// data of hierarchical levels of ever coarser meshes
		class Level
		{
			// number of grid points per dimension
			// at the given moment, only n[0]==n[1] is supported
			int    [] n;
			// reaction coefficient (spatially varying or constant)
			// a = nvar x nvar x naij, where naij is 1xn^2 in case
			// of variable coefficients or 1x1 in case of constant coefficients
			double [][][] a;
			// residual (temporary variable)
			double [][] res;
			// solution
			double [][] x;
			// right hand side
			double [][] b;
			
			public void init(int [] n, int l)
			{
				int nn   = n[0]*n[1];
				this.n   = n;
				this.x   = new double [l][nn];
				this.res = new double [l][nn];
				this.b   = new double [l][nn];
			} // Level::init
			public Level()
			{
			}
		} // class Level
	
		Level [] lvl;

		//
		// alogrithm specific parameters
		//
		// number of smoothing iterations per pre/post-cyle smoothing
		public int nsmooth = 1;
		// number of times the recursive cycle is repeated
		// 1 : v-cycle, 2 : w-cycle
		public int nsubcycle = 1;
		// smoothing factor
		public double o = 2.0/3.0;
		// maximum number of iterations
		public int nmaxiter = 1000;
		// relative tolerance
		public double reltol = 1e-7;

		//
		// problem specific parameters
		//
		// number of (coubled) variables
		int nvar;
		// side length of domain
		double [] L;
		// aux variable indicating which of the nvarxnvar reaction coefficients
		// are scalar (constant)
		private int s_[][];
	
		// diffusion rate (not necessarily isotropic)
		// 2 x nvar
		double [][] d;

		// advection coefficient
		// 2x3xnvar
		double [][][] ad;

		// statistics with infomation about the iteration process
		// numer of iterations
		public int iter;
		// final residual
		public double resn;
		public double relres;
		// covergence flag : 0 converged, -1 max iter reached
		public int flag;

		boolean iszero;

	// default constructor
	public Multigrid_java()
	{
	} // Multigrid

	public void init(double [][][] a, double [][][] ad, double [][] d, double [] L, int [] n) throws Exception
	{
		if (n[0] != n[1])
		{
			throw new Exception("number of grid points in both dimensions must match");
		}

		//this.d = d;
		this.L = L;
		this.nvar = d[0].length;
		this.s_ = new int[this.nvar][this.nvar];
		
		int nl = 1 + (int) Math.round(Math.log(n[0]) / Math.log(2.0));
		if ( (1 << (nl-1)) != n[0])
		{
			throw new Exception("number of grid points must be a power of two");
		}

		// allocate sub-levels
		this.lvl = new Level [nl];
		for (int k=0; k<nl;k++) {this.lvl[k] = new Level();}
		this.lvl[0].n = n;
		this.lvl[0].a = a;
		// only the temporary variable res has to be allocated for the topmost level
		this.lvl[0].res = new double[this.nvar][n[0]*n[1]];
		for (int k=1; k<nl; k++)
		{
			int [] n_ = {this.lvl[k-1].n[0]/2,this.lvl[k-1].n[1]/2};
			this.lvl[k].init(n_, this.nvar);

		} // for k
		// set coefficients
		this.set_coefficients(a,ad,d);
	} // init

	public void set_nmaxiter(int nmaxiter)
	{
		this.nmaxiter = nmaxiter;
	}

	public void set_reltol(double reltol)
	{
		this.reltol = reltol;
	}

	public double [][] get_d()
	{
		return(this.d);
	}

	// set reaction coefficient
	public void set_coefficients(final double [][][] a, final double [][][] ad, final double [][] d)
	{
		this.d  = d;
		this.ad = ad;

		// determine which coefficients are scalars
		for (int i = 0; i < this.nvar; i++)
		{
			for (int j = 0; j < this.nvar; j++)
			{
				if (a[i][j].length > 1)
				{
					this.s_[i][j] = 1;
				} else {
					this.s_[i][j] = 0;
				}
			} // for j
		} // for i
 
		// highest level
		this.lvl[0].a = a;
		for (int k = 1; k<this.lvl.length;k++)
		{
		this.lvl[k].a = new double [this.nvar][this.nvar][];

		for (int i = 0; i < this.nvar; i++)
		{
		for (int j = 0; j < this.nvar; j++)
		{
			if (this.lvl[k-1].a[i][j].length > 1)
			{
				if (this.lvl[k].a[i][j] == null)
				{
					this.lvl[k].a[i][j] = new double[this.lvl[k].n[0]*this.lvl[k].n[1]];
				}
				Resample.down12(this.lvl[k-1].a[i][j],this.lvl[k].a[i][j],this.lvl[k-1].n);
			} else {
				this.lvl[k].a[i][j] = this.lvl[k-1].a[i][j];
			}
		} // j
		} // i
		} // for k
	} // set_a

	public double [][] mldivide(double [][] b, double [][] x)
	{

		int iter = 0;
		this.lvl[0].x = x;
		this.lvl[0].b = b;
		// compute residual res_[0] from x_a[0]
		resfun(0);
		double resn0 = rms(this.lvl[0].res);
		double resn;
		while (true)
		{
			iter = iter+1;
			this.iszero = false;
			cycle(0);
			resfun(0);
			resn = rms(this.lvl[0].res);
			if (resn <= this.reltol*resn0)
			{
				this.flag = 0;
				break;
			}
			if (iter >= this.nmaxiter)
			{
				this.flag = -1;
				break;
			}
			// TODO check for stagnation
		} // while
		this.iter = iter;
		this.resn = resn;
		this.relres = resn/resn0;
		return(this.lvl[0].x);
	} // mldivide

	public double [][] cycle1(final double [][] b)
	{
		if (this.lvl[0].x == null)
		{
			int [] n = this.lvl[0].n;
			this.lvl[0].x = new double [this.nvar][n[0]*n[1]]; 
			//this.lvl[0].res = new double [this.nvar][n[0]*n[1]]; 
		}

		this.lvl[0].b = b;
		this.iszero=true;
		cycle(0);
		return (this.lvl[0].x);
	}

	//public void cycle(double [] b, double [] x,int [] n, int k)
	public void cycle(int k)
	{
		int [] n      = this.lvl[k].n;

		if (n[0]>1)
		{
		// TODO only if nsmooth > 0
		if (this.iszero)
		{
			jacobi_step0(k);
		} else {
			jacobi_step(k);
			this.iszero = true;
		}
		for (int i=1; i<this.nsmooth; i++)
		{
			// smooth, compute x[k] from b, using res[k] as tempvar
			jacobi_step(k);
		} // for i
		for (int i=0; i<this.nsubcycle; i++)
		{
			// residual, compute res_[k] from x_[k] and b_[k]
			resfun(k);
			// downsample, compute b_a[k+1] from res_[k]
			downsample12(this.lvl[k].res,this.lvl[k+1].b,k);
			//downsample12(k);
			// zero initial guess (can go up), and can be skipped
//			for (int l = 0; l<this.nvar;l++)
//			{
//				Arrays.fill(this.lvl[k+1].x[l],0.0);
//				//for (int j = 0;j<n[0]*n[1]/4;j++)
//				//{
//				//	this.lvl[k+1].x[l][j] = 0;
//				//} // for i
//			} // for l
			// cycle, compute x_[k+1] from b_[k+1], using res_[k] as tmp
			cycle(k+1);
			// upsample, compute res_k from x_[k+1]
			upsample12(k);
			// correct x_k from res_k
			for (int l = 0; l<this.nvar; l++)
			{
				for (int j=0; j<n[0]*n[1]; j++)
				{
					//x[j] = x[j] - e[j];
					this.lvl[k].x[l][j] = this.lvl[k].x[l][j] - this.lvl[k].res[l][j]; //e_a[k][j];
				} // for j
			} // for l
		} // for i = 0 .. gamma
		// post-smooth
		for (int i=0; i<this.nsmooth; i++)
		{
			jacobi_step(k);
		} // for i
		} else {
			// final step
			double[][] A = new double[this.nvar][this.nvar];
			for (int i=0; i<this.nvar;i++) 
			{
				for (int j=0; j<this.nvar;j++) 
				{
				A[i][j] = this.lvl[k].a[i][j][0];
				}
			}
			if (this.nvar != 3)
			{
			double[][] b = new double[this.nvar][1];
			for (int i = 0; i<this.nvar; i++)
			{
				b[i][0] = this.lvl[k].b[i][0];
			}
			mldivide_small(A,b,k);
			} else
			{
				double[] b = new double[this.nvar];
				for (int i = 0; i<this.nvar; i++)
				{
					b[i] = this.lvl[k].b[i][0];
				}
				double [] x = Matrix_java.mldivide3(A,b);
				for (int i = 0; i<this.nvar; i++)
				{
					this.lvl[k].x[i][0] = x[i];
				}
			}
		} // else of if n>1
		//return(x);
	} // for cycle

//	public void mldivide_3x3(final int k)
//	{
//		
//	}

	public void mldivide_small(double [][] A, double [][] b,final int k)
	{
		// TODO we can do this right during setup of aa
		Matrix Am = new Matrix(A);
		Matrix bm = new Matrix(b);
		Matrix x = Am.solve(bm);
		for (int i = 0; i<this.nvar; i++)
		{
			this.lvl[k].x[i][0] = x.get(i,0);
		}
		// TODO solve properly nxn
		//for (int l=0; l<this.nvar;l++)
		//{
		//	this.lvl[k].x[l][0] = this.lvl[k].b[l][0]/this.lvl[k].a[l][l][0];
		//}
	}

	//public double [] cycle_(double [] b, double [] x,int [] n, int k)
	//{
	//	return(x);
	//}
	public void downsample12(double [][] in, double [][] out, final int k)
	{
		for (int l=0; l<this.nvar; l++)
		{
			Resample.down12(in[l],out[l],this.lvl[k].n);
		}
		//double [][] tmp = downsample1(in,k);
		//downsample2(tmp,out,k);
		//return(x);
	} // downsample12

	public double [][] downsample12_(double [][] res, int k)
	{
		this.lvl[k].res = res;
		//downsample12(k);
		return(this.lvl[k+1].b);
	} // downsample12

	public void upsample12(int k)
	{
		int [] n = this.lvl[k+1].n;
		for (int l=0; l<this.nvar;l++)
		{
			double [] in = this.lvl[k+1].x[l];
			double [] out = this.lvl[k].res[l];
			Resample.up12(in,out,n);
		}
	} // upsample12

	public double [][] jacobi_step_(double [][] b, double [][] x, int k)
	{
		this.lvl[k].x = x;
		this.lvl[k].b = b;
		jacobi_step(k); //res, b, x, n);
		return(this.lvl[k].x);
	}
	
	//public void jacobi_step(double [] res, double [] b, double [] x, int [] n)
	public void jacobi_step(int k)
	{
		// diagonal of the matrix
		int []  n = this.lvl[k].n;
		double [] dx = {this.L[0]/n[0], this.L[1]/n[1]};
		resfun(k);
		for (int l=0; l<this.nvar; l++)
		{
			if (this.lvl[k].a[l][l].length > 1)
			{
				double d = (   this.d[0][l]*2.0/(dx[0]*dx[0])
                                             + this.d[1][l]*2.0/(dx[1]*dx[1])
					     - this.ad[0][1][l]/dx[0] 
					     - this.ad[1][1][l]/dx[1]
					   );	
				for (int i=0; i< n[0]*n[1];i++)
				{
					this.lvl[k].x[l][i] = this.lvl[k].x[l][i] - this.o*this.lvl[k].res[l][i]/(this.lvl[k].a[l][l][i] + d);
				} // for i
			} else {
				double di   = this.o/(   this.lvl[k].a[l][l][0] 
						       + this.d[0][l]*2.0/(dx[0]*dx[0]) 
						       + this.d[1][l]*2.0/(dx[1]*dx[1])
					     	       - this.ad[0][1][l]/dx[0] 
						       - this.ad[1][1][l]/dx[1]
					             );	
				for (int i=0; i<n[0]*n[1];i++)
				{
					this.lvl[k].x[l][i] = this.lvl[k].x[l][i] - di*this.lvl[k].res[l][i];
				} // for i
			}
		} // for l
	} // jacobi_step

	//public void jacobi_step(double [] res, double [] b, double [] x, int [] n)
	public void jacobi_step0(int k)
	{
		// diagonal of the matrix
		int []  n = this.lvl[k].n;
		double [] dx = {this.L[0]/n[0], this.L[1]/n[1]};
		resfun0(k);
		for (int l=0; l<this.nvar; l++)
		{
			if (this.lvl[k].a[l][l].length > 1)
			{
				double d = (   this.d[0][l]*2.0/(dx[0]*dx[0])
					     + this.d[1][l]*2.0/(dx[1]*dx[1])
					     - this.ad[0][1][l]/dx[0]
					     - this.ad[1][1][l]/dx[1] 
					   );
				for (int i=0; i< n[0]*n[1];i++)
				{
					this.lvl[k].x[l][i] = - this.o*this.lvl[k].res[l][i]/(this.lvl[k].a[l][l][i] + d);
				} // for i
			} else {
				// TODO di can be precomputed
				double di   = this.o/(this.lvl[k].a[l][l][0] 
						       + this.d[0][l]*2.0/(dx[0]*dx[0]) 
						       + this.d[1][l]*2.0/(dx[1]*dx[1])
					     	       - this.ad[0][1][l]/dx[0] 
						       - this.ad[1][1][l]/dx[1]
					             );	
				for (int i=0; i<n[0]*n[1];i++)
				{
					this.lvl[k].x[l][i] = - di*this.lvl[k].res[l][i];
				} // for i
			}
		} // for l
	} // jacobi_step0

	public void resfun0(final int k)
	{
		int [] n      = this.lvl[k].n;
		for (int l=0; l<this.nvar; l++)
		{
			double [] b   = this.lvl[k].b[l];
			double [] res = this.lvl[k].res[l];
			for (int i = 0; i<n[0]*n[1];i++)
			{	
				res[i] = -b[i];
			}
		}
	}

	public double [][] resfun_(double [][] b, double [][] x, final int k)
	{
		this.lvl[k].x = x;
		this.lvl[k].b = b;
		resfun(k);
		return(this.lvl[k].res);
	}

	// note that dimensions "left/right" and "up/down" differ between matlab and java
	public void resfun(int k)
	{
		int [] n      = this.lvl[k].n;
		double [] dx  = {this.L[0]/n[0], this.L[1]/n[1]};
		double xx [][] = this.lvl[k].x;

		for (int l=0; l<this.nvar; l++)
		{
			double [] b   = this.lvl[k].b[l];
			double [] res = this.lvl[k].res[l];
			double [] x   = this.lvl[k].x[l];
			double [][] a = this.lvl[k].a[l];
			
			// off-diagonal elements
			//double r  = (this.d[l]*1.0/(dx[0]*dx[1]));
			// diagonal elements
			// double di  = (4.0*r);
			
			// center
			double cc = (   this.d[0][l]*2.0/(dx[0]*dx[0]) 
				      + this.d[1][l]*2.0/(dx[1]*dx[1])
				      - this.ad[0][1][l]/dx[0]
				      - this.ad[1][1][l]/dx[1]
				    );
			// left
			double cl = ( - this.d[0][l]/(dx[0]*dx[0])
				      - this.ad[0][0][l]/dx[0]
				    );
			// right
			double cr = ( - this.d[0][l]/(dx[0]*dx[0])
				      - this.ad[0][2][l]/dx[0]
				    );
			// up
			double cu = ( - this.d[1][l]/(dx[1]*dx[1])
				      - this.ad[1][0][l]/dx[1]
				    );
			// down
			double cd = ( - this.d[1][l]/(dx[1]*dx[1])
				      - this.ad[1][2][l]/dx[1]
				    );

		// left top corner i=0,j=0
				res[0] = (   cc*x[0]
					   + cl*x[n[1]-1]
					   + cr*x[1]
					   + cu*x[n[0]*n[1]-n[1]]
					   + cd*x[n[1]] 
					   - b[0] 
					 );
				// reaction part
				for (int ll=0; ll<this.nvar;ll++){ res[0] += a[ll][0]*xx[ll][0]; }
		// top right corner i=0,j=n[1]-1
				res[(n[1]-1)] = ( cc*x[n[1]-1]
						+  cl*x[n[1]-2]
						+ cr*x[0]
						+ cu*x[n[0]*n[1]-n[1]+(n[1]-1)]
						+ cd*x[n[1]+(n[1]-1)]
						- b[(n[1]-1)] 
						);
				for (int ll=0; ll<this.nvar;ll++){ res[n[1]-1] += a[ll][s_[l][ll]*((n[1]-1))]*xx[ll][n[1]-1]; }
		

		// first row, i = 0
			for (int j = 1; j<n[1]-1; j++)
			{
				res[j] = (   cc*x[j]
					   + cl*x[j-1]
					   + cr*x[j+1]
				           + cu*x[n[0]*n[1]-n[1]+j]
					   + cd*x[n[1]+j]
					   - b[j] 
					);
				 for (int ll=0; ll<this.nvar;ll++){ res[j] += a[ll][s_[l][ll]*j]*xx[ll][j]; }
			} // for j

		for (int i = 1; i<n[0]-1; i++)
		{
			// first column j == 0
				res[i*n[1]] = (   cc*x[i*n[1]]
						+ cl*x[i*n[1]+n[1]-1]
						+ cr*x[i*n[1]+1]
						+ cu*x[(i-1)*n[1]]
						+ cd*x[(i+1)*n[1]]
						- b[i*n[1]] 
						);
				 for (int ll=0; ll<this.nvar;ll++){ res[i*n[1]] += a[ll][s_[l][ll]*(i*n[1])]*xx[ll][i*n[1]]; }
			for (int j = 1; j<n[1]-1; j++)
			{
				res[i*n[1]+j] = ( cc*x[i*n[1]+j]
						    + cl*x[i*n[1]+j-1]
						    + cr*x[i*n[1]+j+1]
						    + cu*x[(i-1)*n[1]+j]
						    + cd*x[(i+1)*n[1]+j]
						- b[i*n[1]+j] 
						);
				 for (int ll=0; ll<this.nvar;ll++){ res[i*n[1]+j] += a[ll][s_[l][ll]*(i*n[1]+j)]*xx[ll][i*n[1]+j]; }
			} // for j
			// last column, j = n1-1
				res[i*n[1]+(n[1]-1)] = ( cc*x[i*n[1]+n[1]-1]
						    + cl*x[i*n[1]+n[1]-2]
						    + cr*x[i*n[1]+0]
						    + cu*x[(i-1)*n[1]+n[1]-1]
						    + cd*x[(i+1)*n[1]+n[1]-1]
						- b[i*n[1]+(n[1]-1)] 
						);
				 for (int ll=0; ll<this.nvar;ll++){ res[i*n[1]+(n[1]-1)] += a[ll][s_[l][ll]*(i*n[1]+(n[1]-1))]*xx[ll][i*n[1]+(n[1]-1)]; }
		} // for i

		// last column i = n[0]-1
			for (int j = 1; j<n[1]-1; j++)
			{
				res[(n[0]-1)*n[1]+j] = ( cc*x[(n[0]-1)*n[1]+j]
						    + cl*x[(n[0]-1)*n[1]+j-1]
						    + cr*x[(n[0]-1)*n[1]+j+1]
						    + cu*x[(n[0]-2)*n[1]+j]
						    + cd*x[j]
						- b[(n[0]-1)*n[1]+j] 
						);
				 for (int ll=0; ll<this.nvar;ll++){ res[(n[0]-1)*n[1]+j] += a[ll][s_[l][ll]*((n[0]-1)*n[1]+j)]*xx[ll][(n[0]-1)*n[1]+j]; }
			} // for j
		// bottom left corner, i=n[0]-1,j=0
				res[(n[0]-1)*n[1]] = ( cc*x[(n[0]-1)*n[1]]
						    + cl*x[(n[0]-1)*n[1]-1+n[1]]
						    + cr*x[(n[0]-1)*n[1]+1]
						    + cu*x[(n[0]-2)*n[1]]
						    + cd*x[0]
						- b[(n[0]-1)*n[1]]
						); 
				 for (int ll=0; ll<this.nvar;ll++){ res[(n[0]-1)*n[1]] += a[ll][s_[l][ll]*(n[0]-1)]*xx[ll][(n[0]-1)*n[1]]; }
		// bottom left corner, i=n[0]-1,j=n[1]-1
				res[(n[0]-1)*n[1]+(n[1]-1)] = ( cc*x[(n[0]-1)*n[1]+(n[1]-1)]
						    + cl*x[(n[0]-1)*n[1]+(n[1]-2)]
						    + cr*x[(n[0]-1)*n[1]+(0)]
						    + cu*x[(n[0]-2)*n[1]+(n[1]-1)]
						    + cd*x[(n[1]-1)]
						- b[(n[0]-1)*n[1]+(n[1]-1)] 
						);
				 for (int ll=0; ll<this.nvar;ll++){ res[ (n[0]-1)*n[1]+(n[1]-1)]+= a[ll][s_[l][ll]*((n[0]-1)*n[1]+(n[1]-1))]*xx[ll][(n[0]-1)*n[1]+(n[1]-1)]; }
		// return(res);
		} // for vdx
	} // resfun


	public static double rms(final double [][] x)
	{
		double ss = 0;
		int n = 0;
		for (int l=0; l<x.length; l++)
		{
			n = n+x[l].length;
			for (int i=0; i<x[l].length; i++)
			{
				ss = ss + x[l][i]*x[l][i];
			} // for i
		} // for l
		return (Math.sqrt(ss/n));
	} // rms


} // class Multigrid

