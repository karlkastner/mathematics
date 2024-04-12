// 2024-02-01 10:49:41.444719979 +0100
class Resample
{
	public Resample()
	{
	}

	public static double [] down12(double [] in, final int [] n)
	{
		double [] out = new double[n[0]*n[1]/4];
		//for (int i=0; i<10; i++)
		//{
		down12(in, out, n);
		//}
		return(out);
	}
	public static double [] down12_(double [] in, final int [] n)
	{
		double [] out = new double[n[0]*n[1]/4];
		//for (int i=0; i<10; i++)
		//{
		down12_(in, out,n);
		//}
		return (out);
	}

	public static void down12_(double [] in, double [] out, final int [] n)
	{
		double [] tmp = new double [n[0]*n[1]/2];
		Resample.down1(in, tmp, n);
		Resample.down2(tmp, out, n);
	} // resample12

	public static double [] down1(double [] in, final int [] n)
	{
		double [] out = new double[n[0]*n[1]/2];
		down1(in, out, n);
		return(out);
	}
	public static void down12(double [] in, double out [], final int [] n)
	{
		double [][] tmp = new double[3][n[1]/2];
		double [] tmp0 = tmp[0];
		// first row

		// last row, first element
		tmp[0][0] = (     0.50*in[(n[0]-1)*n[1]]
		              + 0.25*( in[(n[0]-1)*n[1]+n[1]-1] 
		             +         in[(n[0]-1)*n[1]+1]
			     )
			   );

		// last row, i = 0
		for (int j = 1; j<n[1]/2; j++)
		{
			tmp[0][j] = (     0.50*in[(n[0]-1)*n[1]+2*j]
			              + 0.25*( in[(n[0]-1)*n[1]+2*j-1] 
			             +         in[(n[0]-1)*n[1]+2*j+1]
				     )
				   );
		} // for j
		// remaining rows
		for (int i = 0; i<n[0]/2; i++)
		{
			// first column, j = 0
			tmp[2][0] = (   0.50*in[(2*i+1)*n[1] + 0]
			            + 0.25*( in[(2*i+1)*n[1] + n[1] - 1] 
			                   + in[(2*i+1)*n[1] + 1]
						     )
					);
			tmp[1][0] = (   0.50*in[2*i*n[1] + 0]
			            + 0.25*( in[2*i*n[1] + n[1] - 1] 
			                   + in[2*i*n[1] + 1]
						     )
					);
			out[(i)*n[1]/2 + 0] = (      0.50*tmp[1][0]
			                      + 0.25*(  tmp[0][0] 
			                             +  tmp[2][0]
						     )
					    );

			for (int j = 1; j<n[1]/2; j++)
			{	
			// first dimension (down2)
			tmp[2][j] = (     0.50*in[(2*i+1)*n[1] + 2*j]
			              + 0.25*( in[(2*i+1)*n[1] + 2*j-1] 
			             +         in[(2*i+1)*n[1] + 2*j+1]
				     )
				   );
			tmp[1][j] = (     0.50*in[2*i*n[1] + 2*j]
			              + 0.25*( in[2*i*n[1] + 2*j-1] 
			             +         in[2*i*n[1] + 2*j+1]
				     )
				   );
			// second dimension
			out[i*n[1]/2 + j] = (      0.50*tmp[1][j]
			                             + 0.25*(  tmp[0][j] 
			                                    +  tmp[2][j]
						     )
					    );
			} // for j
			// circle
			double [] aux = tmp[0];
			tmp[0] = tmp[2];
			//tmp_[1] = tmp[2];
			tmp[2] = aux;

		} // for i
	} // down12_

	public static void down1(double [] in, double out [], final int [] n)
	{
		// downsample along first dimension
		for (int j = 0; j<n[1]; j++)
		{
		// first element, i == 0
			out[j] = (   0.50*in[j]
			                      + 0.25*( in[(n[0]-1)*n[1] + j] 
			                             + in[ n[1] + j]
						     )
					    );
		}

		for (int i = 1; i<n[0]/2; i++)
		{
		for (int j = 0; j<n[1]; j++)
		{
			out[i*n[1] + j] = (      0.50*in[2*i*n[1] + j]
			                      + 0.25*( in[(2*i-1)*n[1] + j] 
			                             + in[(2*i+1)*n[1] + j]
						     )
					    );
		} // for j
		} // for i
	}

	public static void down2(double [] x, double [] x_, final int [] n)
	{
		// first column
		for (int i = 0; i<n[0]/2; i++)
		{
			x_[i*n[1]/2 + 0] = (   0.50*x[i*n[1] + 0]
			                      + 0.25*( x[i*n[1] + n[1] - 1] 
			                             + x[i*n[1] + 1]
						     )
					);
		}
	
		// downsample along first dimension
		for (int i = 0; i<n[0]/2; i++)
		{
		for (int j = 1; j<n[1]/2; j++)
		{
			x_[i*n[1]/2 + j] = (   0.50*x[i*n[1] + 2*j]
			                      + 0.25*( x[i*n[1] + 2*j-1] 
			                             + x[i*n[1] + 2*j+1]
						     )
					    );
		} // for j
		} // for i
	} // down2

	public static double [] up1(double [] x, int [] n) 
	{
		double [] tmp = new double[2*n[0]*n[1]];

		for (int i = 0; i<n[0]-1; i++)
		{
		for (int j = 0; j<n[1]; j++)
		{
			// even, take over
			tmp[2*i*n[1] + j]     = x[i*n[1] + j];
			// odd, interpolate
			tmp[(2*i+1)*n[1] + j] = 0.5*(x[i*n[1] + j] + x[(i+1)*n[1] + j]);
		}}

		for (int j = 0; j<n[1]; j++)
		{	
		// last column
			tmp[2*(n[0]-1)*n[1] + j]     = x[(n[0]-1)*n[1] + j];
			// odd, interpolate
			tmp[(2*n[0]-1)*n[1] + j] = 0.5*(x[(n[0]-1)*n[1] + j] + x[0 + j]);
		}
		return(tmp);
	} // up1

	public static void up2(final double [] tmp, final double [] out, final int [] n) //k)
	{
		for (int i = 0; i<2*n[0]; i++)
		{
		for (int j = 0; j<n[1]-1; j++)
		{
			out[i*(2*n[1]) + 2*j] = tmp[i*n[1] + j]; 
			out[i*(2*n[1]) + 2*j+1] = 0.5*(tmp[i*n[1] + j] + tmp[i*n[1] + j + 1]); 
		} // for j
		// last row
			out[i*(2*n[1]) + 2*(n[1]-1)] = tmp[i*n[1] + n[1]-1]; 
			out[i*(2*n[1]) + 2*(n[1]-1)+1] = 0.5*(tmp[i*n[1] + n[1]-1] + tmp[i*n[1] + 0]); 
		} // for i
	}

	public static void up12(double [] in, double [] out, int [] n)
	{
		// top left i == 0, j == 0
		out[4*0*n[1]+ 2*0] = in[0*n[1]+0];
		// top right
		out[4*0*n[1]+ 2*n[1]-1] = 0.5*(in[0*n[1]+ 0] +in[0*n[1]+ n[1]-1]);
		// bottom left
		out[(4*n[0]-2)*n[1]+ 2*0] = 0.5*(in[(n[0]-1)*n[1]+ 0] +in[0*n[1]+ 0]);


		// i == 0, first column
		for (int j=1;j<n[1];j++)
		{
			// set value
			out[4*0*n[1]+ 2*j] = in[0*n[1]+j];
			// interpolate in row
			out[4*0*n[1]+ 2*j-1] = 0.5*(in[0*n[1]+ j] +in[0*n[1]+ j-1]);
		}


		for (int i=1; i<n[0];i++)
		{
			// j == 0, first row
			out[4*i*n[1]+ 2*0] = in[i*n[1]+0];
			out[4*i*n[1]+ 2*n[1]-1] = 0.5*(in[i*n[1]+ 0] +in[i*n[1]+ n[1]-1]);
			out[(4*i-2)*n[1]+ 2*0] = 0.5*(in[(i-1)*n[1]+ 0] +in[i*n[1]+ 0]);
			out[(4*i-2)*n[1]+ 2*n[1]-1] = 0.5*(out[4*i*n[1] + 2*n[1]-1] + out[(4*i-4)*n[1]+2*n[1]-1]);
			
		// inner
		for (int j=1;j<n[1];j++)
		{
			// set value
			out[4*i*n[1]+ 2*j] = in[i*n[1]+j];
			// interpolate in row
			out[4*i*n[1]+ 2*j-1] = 0.5*(in[i*n[1]+ j] +in[i*n[1]+ j-1]);
			// intepolate column
			out[(4*i-2)*n[1]+ 2*j] = 0.5*(in[(i-1)*n[1]+ j] +in[i*n[1]+ j]);
			// interpolate cross (note that there are two alternative ways to cross-interpolate)
			out[(4*i-2)*n[1]+ 2*j-1] = 0.5*(out[4*i*n[1] + 2*j-1] + out[(4*i-4)*n[1]+2*j-1]);
			//0.5*(in[i-1, j] +in[i, j]);
		}
		}
		// last column
		for (int j=1;j<n[1];j++)
		{
			// intepolate last column
			out[(4*n[0]-2)*n[1]+ 2*j] = 0.5*(in[(n[0]-1)*n[1]+ j] + in[0*n[1]+ j]);
			// interpolate cross (note that there are two alternative ways to cross-interpolate)
			out[(4*n[0]-2)*n[1]+ 2*j-1] = 0.5*(out[4*0*n[1] + 2*j-1] + out[(4*n[0]-4)*n[1]+2*j-1]);
			//0.5*(in[i-1, j] +in[i, j]);
		}
		// bottom right corner
		out[(4*n[0]-2)*n[1]+ 2*n[1]-1] = 0.5*(out[4*0*n[1] + 2*n[1]-1] + out[(4*n[0]-4)*n[1]+2*n[1]-1]);

	} // up12

	public static void up12_(double [] in, double [] out, int [] n)
	{
		double [] tmp = up1(in,n);
		up2(tmp,out,n);
	}

	public static double [] up12_(double [] in, int [] n)
	{
		double [] out = new double [4*n[0]*n[1]];
		//for (int i = 0; i<1000; i++)
		//{
		double [] tmp = up1(in,n);
		up2(tmp,out,n);
		//}
		return(out);
	}


	public static double [] up12(double [] in, int [] n)
	{
		double [] out = new double [4*n[0]*n[1]];
		//for (int i = 0; i<1000; i++)
		//{
			up12(in,out,n);
		//}
		return(out);
	}

//	public double [][] up12_(double [][] x, int k)
	//{
	//	this.lvl[k+1].x = x;
	//	up12(k);
	//	return(this.lvl[k].res);
	//}



} // resample
