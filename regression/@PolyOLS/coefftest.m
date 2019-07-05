% 2015-06-08 19:32:53.298194708 +0200
function obj = coefftest(obj)
	% note : columns should be orthogonal
	nsigma = 2;
	c      = tinv(normcdf(nsigma),obj.nsample-obj.nparam);
	d      = c*diag(obj.C)
	[obj.param - d, obj.param, obj.param+d, abs(obj.param)>abs(d), d./abs(obj.param) ]
end
