% 2015-02-19 17:47:55.246841347 +0100
% Karl Kastner, Berlin
%
% function w = hanwin(x,x0,range)
%
%% hanning filter window
% effective sample size 2/3 n for uncorrelated data
%
% function w = hanwin(x,x0,L)
function w = hanwin(x,x0,L)
	if (nargin() < 2 || isempty(x0))
		x0 = 0.5*(x(end)+x(1));
	end
	x = x-x0;
	if (nargin() < 3 || isempty(L))
		L = x(end);
		% scale up such that the end values or non zero
		n = length(x);
		L = L*n/(n-1);
	end
	x = x/(2*L);
	w = 0.5*(1+cos(2*pi*x)).*(x >= -0.5).*(x <= 0.5);
	w = w/sum(w);
end

