% Mo 15. Feb 16:54:19 CET 2016
% Karl Kastner, Berlin
% if frequency is higher than 
function [rho, L, neff] = cosexpdelay(acf,n,ds)
	if (nargin() < 2)
		n = 1;
	end
	if (nargin() < 3)
		ds = 1;
	end
	[rho L] = ar1delay(acf);
	%p0 = [1 1];
	p0 = real([1./L,1./L]);

	x = ds*(0:length(acf)-1)';
	a = cvec(acf);

	fun = @(p) cos(2*pi*x*p(1)).*exp(-x*p(2));
	opt = struct();
	opt.Display = 'off';
	p   = lsqnonlin(@(p) (fun(p)-acf),p0,[],[],opt);
	% TODO, do this analyticaly
	f = fun(p);
	L = 0.5*ds*(f(1) + 2*sum(f(2:end)));
	neff = (ds*n)/(2*L);
%	L
%	1./p
%pause
%	plot([a f])
%pause

%	L = 1./p(2);
	rho = exp(-1/L);

	


