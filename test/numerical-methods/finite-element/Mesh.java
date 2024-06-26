// Sat Mar 15 19:30:16 WIB 2014
// Karl Kastner, Berlin

// super class for meshes in n-dimensions
abstract class Mesh
{
	// point coordinates
	// points are not only corner vertices, but also other points of an element
	// for higher order basis functions
	// TODO either split this up to X Y Z or make first dimension the smaller one
	public double [][] P;
	// number of points (nt >= P.length)
	public int np;

	// elements (line segments, triangles, quads, tetras, hexas, ...)
	// defined by indices into rows of P
	public int [][] T;
	// number of triangles (nt >= T.length)
	public int nt;

//	// alternative triagulation with basis function p/2 for error estimation
//	public int [][] S;
//	public int ns;

	// boundary elments, defined by indices into rows of P
	public int [][] Bc;
	// number of boundary elements (nb >= Bc.length)
	public int nb;

	// indices of neighbouring elements
	public int [][] N;

	// dimension
	int DIM; // final

	// tests, whether the point p is contained in the elment with index tdx
	abstract public boolean contains(final int tdx, final double [] p);
} // interface Mesh

