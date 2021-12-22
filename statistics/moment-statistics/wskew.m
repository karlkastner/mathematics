% Do 11. Feb 19:33:24 CET 2016
% Karl Kastner, Berlin
%
%% skewness of a weighted set of samples
% function sk = wskew(w,x)
function sk = wskew(w,x)
	w  = w./sum(w);
	mu = wmean(w,x);
	sd = wstd(w,x);
	c3 = w'*(x-mu).^3;
	sk = c3/(sd^3);
end
