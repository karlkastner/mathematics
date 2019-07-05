% 2016-03-04 11:16:28.144823961 +0100
% Karl Kastner, Berlin
%
%% range of values within a certain range of indices (window)
%
function R = range_window(x,n)
	x0 = x;
	[x, sdx] = sort(x);
	n = max(1,round(n));
	%R = zeros(size(x));
	m = length(x);

	idc = 1:length(x);
	idl = max(1,idc-n);
	idr = min(m,idc+n);
	rr  = x(idr)-x(idc); 
	rl  = x(idc)-x(idl);
	R   = max([rvec(rr); rvec(rl)])';
end

