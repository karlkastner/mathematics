% Do 11. Feb 19:33:38 CET 2016
% Karl Kastner, Berlin
%
%% kurtosis with weighted samples
function ku = wkurt(w,x)
	w  = w/sum(w);
	mu = wmean(w,x);
	sd = wstd(w,x);
	c4 = w'*(x-mu).^4;
	ku = c4/sd^4;
end

