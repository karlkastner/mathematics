% Fri 17 May 10:47:26 CEST 2024
function par = generalized_normcdf_fit(x,c,par)
	x = double(x);
	c = double(c);
	x = cvec(x);
	c = cvec(c);
	if (nargin()<3)
		par = [0,1,1];
	end
	lb = [-inf,0,0];
	par = lsqnonlin(@(p) c - generalized_normcdf(x,p(1),p(2),p(3)),par,lb);
end

