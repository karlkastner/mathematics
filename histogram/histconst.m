% 2016-03-04 12:06:38.510001041 +0100
% Karl Kastner, Berlin

% constant bin size along x
function out = histconst(x,m)
	n = length(x);
	% grid (container boundaries)
	x = sort(x);
	out.edge.x      = x(1) + (x(end)-x(1))*(0:m)/m;
	out.centre.x    = 0.5*(out.edge.x(1:m)+out.edge.x(2:m+1));
	out.cell.length	= diff(out.edge.x);
	% cdf
	sample.cdf = (1:n)/n;
	% cdf linearly varying between edges
	out.edge.cdf   = interp1(x,sample.cdf,out.edge.x,'linear');
	% cdf interpolated to cell centres
	out.centre.cdf = 0.5*(out.edge.cdf(1:m)+out.edge.cdf(2:m+1));
	% pdf
	out.cell.h     = diff(out.edge.cdf);
	% number of samples
	out.cell.count = out.cell.h.*m;
	% pdf constant over the bin
	out.centre.pdf = out.cell.h ./ out.cell.length;

	% compute the l/h shifted histogram
%	out_.edge.x = 

end

