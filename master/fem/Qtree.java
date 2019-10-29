import java.awt.*;
import javax.swing.*;
import java.awt.geom.*;

// build a tree of bubbles containing the elements within the leaves
// TODO this is not yet prepared for elements with curved boundaries
// it assumes no overlapping elements
// TODO add some statistical functions (depth,number of elements, etc)
class Qtree
{
	private final int nvert;
	// index into the elements
	public int [] index;
	// root of the tree
	private Q q;
	private final Mesh mesh;

	// constructor
	public Qtree(Mesh mesh)
	{
		this.mesh = mesh;
		// TODO dangerous if T not initialised
		nvert = mesh.T[0].length;
		index = new int[mesh.nt];
		// prepare the index
		for (int i=0; i<mesh.nt; i++)
		{
			index[i] = i;
		}
		q = new Q(0,mesh.nt-1);
//		for (int i=0; i<mesh.nt; i++)
//		{
//		System.out.println(index[i]);
//		}
	} // Qtree constructor

	public int contains(double [] p)
	{
		return q.contains(p);
	} // contains

	public void plot()
	{
		JFrame frame = new JFrame();
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		//JPanel panel = new MyPanel();
		JComponent panel = new MyPanel();
		frame.add(panel);
		frame.pack();
		//setMinimumSize(getSize());
		Dimension dim = new Dimension(1000, 500);
		frame.setMinimumSize(dim);
		frame.setVisible(true);
	} // Qtree::plot
 
	private class Q
	{
		final int l;
		final int r;
		final double origin [];
		double radius;
		// children
		// TODO, this could also be organised into an array
		Q child1;
		Q child2;
		public Q(final int l, final int r)
		{
			this.l = l;
			this.r = r;
			double [] max = new double[mesh.DIM];
			double [] min = new double[mesh.DIM];
			origin = new double[mesh.DIM];
			for (int jdx=0; jdx<mesh.DIM; jdx++)
			{
				// TODO sloppy inf
				max[jdx] = -1.0/0.0;
				min[jdx] = +1.0/0.0;
			}

			// get the centre coordinate and radius of the current bubble
			// for each contained element
			for (int i=l; i<=r; i++)
			{
				// for each corner vertex of the element
				for (int kdx=0; kdx<nvert;kdx++)
				{
					// for each dimension
					for (int jdx=0; jdx<mesh.DIM; jdx++)
					{
						max[jdx] = Math.max(max[jdx], mesh.P[mesh.T[index[i]][kdx]-1][jdx]);
						min[jdx] = Math.min(min[jdx], mesh.P[mesh.T[index[i]][kdx]-1][jdx]);
					} // for jdx
				} // for kdx
			} // for idx
			origin[0] = 0.5*(max[0]+min[0]);
			double r2    = (origin[0]-min[0])*(origin[0]-min[0]);
			for (int i=1; i<mesh.DIM; i++)
			{
				origin[i] = 0.5*(max[i]+min[i]);
				//radius = Math.max(radius, origin[i]-min[i]);
				r2 += (origin[i]-min[i])*(origin[i]-min[i]);
			} // for i
			radius = Math.sqrt(r2);
			//System.out.println(origin[0] + " " + origin[1] + " " + radius );

			// split in two if there is more than one contained object
			if (l<r)
			{
				int ll = l;
				int rr = r;
				int sdim = 0;
				// choose widest axis
				for (int i=1; i<mesh.DIM; i++)
				{
					if (origin[i]-min[i] >= origin[sdim]-min[sdim])
					{
						sdim = i;
					} // if
				} // for i
				// distribute contained objects among children
				int lastmove = 0;
				while (ll < rr)
				{
					// get distance from centre of gravity along major axis
					double d = mesh.P[mesh.T[index[ll]][0]-1][sdim];
					for (int i=1; i<nvert; i++)
					{
						d+=mesh.P[mesh.T[index[ll]][i]-1][sdim];
					} // for i
					
					if ( 	d <= nvert*origin[sdim] )
					{
						// object belongs to first child
						ll++;
						lastmove = 0;
					} else {
						// object belongs to right child
						Swap.swap(index,ll,rr);
						rr--;
						lastmove = 1;
					} // if
				} // while
				//System.out.println(l + " " + (ll-1+lastmove) + " " + (rr+lastmove) + " " + r);
				child1 = new Q(l,ll-1+lastmove);
				child2 = new Q(rr+lastmove,r);
			} else {
				child1 = null;
				child2 = null;
			} // if l<r
		} // Q constructor

		public int contains(final double [] p)
		{
			int retval = -1;
			// check if the point is inside the current bubble
			double d2 = (p[0]-origin[0])*(p[0]-origin[0]);
			for (int i=1; i<mesh.DIM; i++)
			{
				d2 += (p[i]-origin[i])*(p[i]-origin[i]);
			} // for int i

			if ( d2 <= radius*radius)
			{
				// check if this is a leaf
				if (l==r)
				{
					// check if the contained object also contains the point
//					System.out.println(l);
//					System.out.println(index[l]);
					if (mesh.contains(index[l],p))
					{
						retval = index[l];
					}
				} else {
					// this is not a leaf, recurse
					retval = child1.contains(p);
					if (retval < 0)
					{
						retval = child2.contains(p);
					}
				}
			} // if bubble contains point
			return retval;
		} // contains

		public void plot(Graphics2D g2, Dimension dim, double [] origin, double radius)
		{
			double [] o = { origin[0]-radius, origin[1]-radius };
			double [] s = { 0.5*dim.width/radius,
					0.5*dim.height/radius };
			g2.draw(new Ellipse2D.Double(
				 (((this.origin[0]-this.radius) - o[0]) * s[0]),
				 (((this.origin[1]-this.radius) - o[1]) * s[1]),
				 (this.radius*s[0]),
				 (this.radius*s[1])));
			g2.drawRect(
				(int) (((this.origin[0]-this.radius) - o[0]) * s[0]),
				(int) (((this.origin[1]-this.radius) - o[1]) * s[1]),
				(int) (this.radius*s[0]),
				(int) (this.radius*s[1]));
//			System.out.println(
//				((origin[0]-radius) - 0) / 1 + " " +
//				((origin[1]-radius) - 0) / 1 + " " +
//				radius / 1 + " " +
//				radius / 1);
			
			if (null != child1)
			{
				child1.plot(g2, dim, origin, radius);
			}
			if (null != child2)
			{
				child2.plot(g2, dim, origin, radius);
			}
		} // plot()

	} // class Q

	public class MyPanel extends JComponent
	{
		@Override
		public void paintComponent(Graphics g)
		{
	    		super.paintComponent(g);
		        Dimension dim = getSize();
		        g.setColor(Color.red);
		        g.setColor(Color.black);
//			System.out.println("Starting recursion");
			q.plot((Graphics2D) g, dim, q.origin, q.radius);
		}
	} // class MyPanel

} // class Qtree

