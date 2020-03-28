% 2016-03-04 12:08:27.950705236 +0100
% Karl Kastner, Berlin
% adaptive bin size in x (constant intervals along p)
function out = histadapt(x,m)
	n = length(x);
	% sort
	x = sort(x);
	% create cdf
	sample.cdf = (1:n)/(n+1);
	% 
	out.edge.cdf  = (0:m)/m;
	out.edge.x    = interp1(sample.cdf,x,out.edge.cdf,'linear','extrap');
	out.centre.x  = 0.5*(out.edge.x(1:m)+out.edge.x(2:m+1));
	out.centre.x = [2*out.centre.x(1)-out.centre.x(2),out.centre.x,2*out.centre.x(end)-out.centre.x(end-1)];

	out.cell.count  = n/m;
	out.cell.h      = 1/m;
	out.cell.length = (diff(out.edge.x));
	out.centre.pdf  = out.cell.h./out.cell.length;
	out.centre.pdf  =[0,out.centre.pdf,0];
end

