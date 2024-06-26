% 2024-05-17 08:56:14.963090742 +0200
% Karl Kastner, Berlin
function par = fit_rational_sigmoid_cdf(x,c,par)
	x = double(x);
	c = double(c);
	x = cvec(x);
	c = cvec(c);
	if (nargin()<3)
		par = [0,1,1];
	end
	par = lsqnonlin(@(p) c - rational_sigmoid_cdf(x,p(1),p(2),p(3)),par);
end
