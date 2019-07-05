// Wed Jun 25 14:04:35 WIB 2014
// Karl Kastner, Berlin
// convert cartesian to streamwise coordinates and back

import java.lang.*;

class XY2sn
{
	// variables for data exchange with matlab
	public double [] S;
	public double [] N;
	public int [][] Nkey;

	// default constructor
	public XY2sn() {}

	// final causes also problems in matlab
	public void xy2sn_(final double [] cX,
			   final double [] cY,
		           final double [] cS,
			   final int    [] cSeg,
                           final double [] X,
                           final double [] Y,
                           double [] S,
                           double [] N) throws Exception
	{
//System.out.println(array.getClass());
		// aux variables
		// find nearest neighbour
		Qtree2 qtree      = new Qtree2(cX,cY);
		int [][] Nkey    = new int[X.length][1];
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
			if (-1 == mx1)
			{
				// no neighbour found likely reason is a NaN coordinate of the current point
				S[idx] = Double.NaN;
				N[idx] = Double.NaN;
				continue;
			} else if (0 == mx1 || cSeg[mx1-1] != cSeg[mx1])
			{
				// first element or first element in this segment
				if (cX.length-1 != mx1 && cSeg[mx1] == cSeg[mx1+1])
				{
					// mx1 is leftmost element,
					// and there is a right neighbour, choose it
					mx2 = mx1+1;
				} else {
					// this segment is of unit length, 
					// projection impossible
					S[idx] = Double.NaN;
					N[idx] = Double.NaN;
					continue;
				}
			} else if(cX.length-1 == mx1 || cSeg[mx1] != cSeg[mx1+1])
			{
				// last element or last element in this segment
				if (0 != mx1 && cSeg[mx1-1] == cSeg[mx1])
				{	
					// mx1 is rightmost element,
					// and there is a left neighbour, choose it
					mx2 = mx1-1;
				} else {
					// this segment is of unit length, 
					// projection impossible
					S[idx] = Double.NaN;
					N[idx] = Double.NaN;
					continue;
				}
			} else {
			// TODO this should not be the second closest but the second neighbour leading to a convex S
				double d2left  = dist2(X[idx],Y[idx],cX[mx1-1],cY[mx1-1]);
				double d2right = dist2(X[idx],Y[idx],cX[mx1+1],cY[mx1+1]);

				if (d2left < d2right)
				{
					mx2 = mx1-1;
				} else {
					mx2 = mx1+1;
				}
			}
			// corner vertices of the triangle
			double Ax = cX[mx1];
			double Ay = cY[mx1];
			double Bx = cX[mx2];
			double By = cY[mx2];
			double Cx = X[idx];
			double Cy = Y[idx];
			// translate 
			// make A the local coodinate origin
			Bx = Bx-Ax;
			By = By-Ay;
			Cx = Cx-Ax;
			Cy = Cy-Ay;

			// rotate
			double hyp = Math.hypot(Bx,By); 
			double sina = (1.0/hyp)*By;
			double cosa = (1.0/hyp)*Bx;
			double Cx_ =  cosa*Cx + sina*Cy;
			double Cy_ = -sina*Cx + cosa*Cy;

			if (cS[mx1] < cS[mx2])
			{
				S[idx] = cS[mx1]+Cx_;
				N[idx] = Cy_;
			} else {
				S[idx] = cS[mx1]-Cx_;
				N[idx] = -Cy_;
			}
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

