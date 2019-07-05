% 2014-12-11 18:49:36.289389338 +0100
% Karl Kastner, Berlin
%
%% regression using the pth-fraction of samples with smallest residual
%
% TODO recursive
% TODO this is not good, determine bad samples by hat matrix
function [c, res2] = areg(x,y)
	p = 0.5;
	n = length(x);
	A = [ones(size(x), x];
	c = A \ y;
	% this should actually be an exterior estimate
	res = A*c - y;
	[void, sdx] = sort(res.*res);
	% only use samples with small residuals
	sdx = sdx(1:ceil(p*n));
	c = A(sdx,:) \ y(sdx);
	res = A*c - y;
	res2 = res'*res/(n-2);	
end

