% Thu 28 Jul 10:45:26 CEST 2016
% Karl Kastner, Berlin
%% normalize columns in x to [-1,1]
function [x,x0,L] = normalize1(x)
	a = min(x);
	b = max(x);
	xc = 0.5*(a+b);
	L  = b-a;
	% x = (x-x0)/L);
	% translate
	x = bsxfun(@minus,x,xc);
	% shift
	x = bsxfun(@times,x,1./L);
end

