% 2016-08-17 16:16:51.155104619 +0200
% Karl Kastner, Berlin
%
%% coefficient of determination when samples are invalid
function r2 = nanr2(x,y)
	f = isfinite(x) & isfinite(y);
	x(~f) = 0;
	y(~f) = 0;
	%if (isempty(fdx))
	%	r2 = NaN;
	%else
		%s2  = mean((x(fdx)-mean(x(fdx))).^2);
		%se2 = mean((y(fdx)-x(fdx)).^2);
		%r2 = 1 - se2/s2;
		n    = sum(f);
		mu   = sum(x)./n;
		mu   = bsxfun(@times,mu,f);
		ss2  = sum((x-mu).^2);
		%sse2 = sum(((y-x).*f).^2);
		%ss2  = sum((x-mu*f).^2);
		sse2 = sum(((y-x)).^2);
		r2 = 1 - sse2./ss2;
	%end
end % r2

