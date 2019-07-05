% 2016-02-15 12:14:21.310572496 +0100
%
%% effective sample size from acf
function [neff T a fdx] = acf_sample_size(a,n,y,dt)
	if (isempty(a))
		n = length(y);
		a = autocorr(y,n-1);
	end
	% 5% change to select a random bin, if there is no correlation
	p = 0.5;

	% probability to reject
	m = length(a);
	x = 1:m;
	P = p./(m-(m-x)*p);
	% z-value for rejection
	Z = norminv(1-0.5*P')/sqrt(n);
%	note: variance of coefficients is not equal!!!
%	Z = 1.96/sqrt(n);
	% this gives some jitter for oscillating acfs, but is acceptable
	%L = find(abs(a(2:end)) > lim)+1;
	fdx = find([1;abs(cvec(a)) > cvec(Z)],1,'last')-1;
	%find(abs(a) > lim,1,'last');
	T = 0.5*a(1)+sum(a(1:fdx));
	neff = n/(2*T);
%	clf
%	plot([a cumsum(a)])
%	vline(fdx)
%	T, pause
	if (nargin() > 3)
		T = dt*T;
	end
end

