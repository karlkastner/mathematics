// Wed Jun 25 14:04:35 WIB 2014
// Karl Kastner, Berlin
// convert cartesian to streamwise coordinates and back
import java.lang.*;

class XY2sn2
{
	// variables for data exchange with matlab
	public double [] S;
	public double [] N;
	public int [][] Nkey;

	// default constructor
	public XY2sn2() {}

	public void xy2sn_(final double [] cX,
		 final double [] cY,
		 final double [] cS,
		 final double [] X,
		 final double [] Y,
		 double [] S,
		 double [] N) throws Exception
	{
		// aux variables
		// find nearest neighbour
		Qtree2 qtree = new Qtree2(cX,cY);
		int [][] Nkey = new int[X.length][1];
		double [][] Dmin = new double[X.length][1]; 
		qtree.nearest_neighbour(X, Y, Nkey, Dmin);

		// for each point in X, find closest centreline point and distance to it
		for (int idx=0; idx<X.length; idx++)
		{
			// get closest centre line point
			int mx1 = Nkey[idx][0]-1;
			int mx2;
			// determine second smallest element
			// take left neighbour if it exists
			// no right neighbour exists
			// or distance to left neighbour is smaller
			// TODO, this may get out of boundx
			double d2left  = dist2(X[idx],Y[idx],cX[mx1-1],cY[mx1-1]);
			double d2right = dist2(X[idx],Y[idx],cX[mx1+1],cY[mx1+1]);

			if (d2left < d2right)
//			if ( (mx1-1>=0)
//			     &&	
//			     ( (mx1+1>=Dmin.length)
//                             || (Dmin[mx1-1][0]<Dmin[mx1+1][0])))
			{
				mx2 = mx1-1;
			} else {
				mx2 = mx1+1;
			}
			// corner vertices of the triangle
			double Ax = cX[mx1];
			double Ay = cY[mx1];
			double Bx = cX[mx2];
			double By = cY[mx2];
			double Cx = X[idx];
			double Cy = Y[idx];
			// make A the local coodinate origin
			Bx = Bx-Ax;
			By = By-Ay;
			Cx = Cx-Ax;
			Cy = Cy-Ay;
			Ax = 0;
			Ay = 0;

			// sides
			double ABx = Bx-Ax;
			double ABy = By-Ay;
			double ACx = Cx-Ax;
			double ACy = Cy-Ay;
			double BCx = Cx-Bx;
			double BCy = Cy-By;
			// squared side length of triangle
			double aa = dot(BCx,BCy,BCx,BCy);
			double bb = dot(ACx,ACy,ACx,ACy);
			double cc = dot(ABx,ABy,ABx,ABy);
			// height of triangle
			// also not optimal take denominator and nominator into one sqrt
			double n = Math.sqrt(2*(aa*bb + bb*cc + aa*cc) - (aa*aa + bb*bb + cc*cc))/(2*Math.sqrt(cc));

			// distance to foot point from mx1
			// ds = bb/sqrt(cc) applies only true to right angled triangles
			// and misses direction, so min ||C-P||, s.t. P = A + alpha AB yields:
			//ds = (C'*AB - A'*AB)/sqrt(cc);
			//P = A+ds*AB/sqrt(cc);
			// s.t. P = pA + (1-p)B
			double p = (dot(Bx,By,Bx,By) - dot(Ax,Ay,Bx,By) + dot(Ax,Ay,Cx,Cy) - dot(Bx,By,Cx,Cy))
				/ (dot(Ax,Ay,Ax,Ay) - 2*dot(Ax,Ay,Bx,By) + dot(Bx,By,Bx,By) );
			double Px = p*Ax + (1-p)*Bx;
			double Py = p*Ay + (1-p)*By;
			double PCx = Cx-Px;
			double PCy = Cy-Py;
			double s = p*cS[mx1] + (1-p)*cS[mx2];
			// n is unsigned, so determine, whether it is left or right of AB
			//n = Math.signum(PCy*ABx*(cS[mx2]-cS[mx1]))*n;
			if (PCy*ABx*(cS[mx2]-cS[mx1]) < 0)
			{
				n = -Math.abs(n);
			}
			S[idx] = s;
			N[idx] = n;
	} // for idx
			// store values for matlab communication
			this.S = S;
			this.N = N;
			this.Nkey = Nkey;
	} // function xy2sn

	private final static double dist2( final double ax, final double ay, final double bx, final double by)
	{
		return (bx-ax)*(bx-ax) + (by-ay)*(by-ay);
	}
			
	private final static double dot(final double ax, final double ay, final double bx, final double by)	
	{
		return (ax*bx)+(ay*by);
	} // function dot
} // class xy2sn

