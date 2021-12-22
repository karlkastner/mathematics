% Wed 28 Jul 17:48:49 CEST 2021
function [s,q] = periodogram_std(fx,S)
	%q = periodogram_quantiles(fx,S,[0.16,0.84],true);
	p = [0.16,0.84];
%	p = [0.25,0.75];
%	p = [0.4,0.6];
	q = periodogram_quantiles(fx,S,p,true);
	r = norminv(p(2))-norminv(p(1));
	%s = sqrt(q(2)/q(1));
	s = (q(2)/q(1)).^(1/r); 
end
