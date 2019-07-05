// 2014-07-10 09:13:28.630762443 +0200
// TODO, this can also be done in matlab as a single bin accumarray
// followed by reshaping
class Bin
{
	public double [][][] bin3_;	

	public Bin() {};

	public void bin3(int [] i1, int [] i2, int [] i3, double [] val, double [][][] bval)
	{
		for (int i=0; i<val.length;i++)
		{
			bval[i1[i]-1][i2[i]-1][i3[i]-1] += val[i];
		}
		bin3_ = bval;
	}
	
}

