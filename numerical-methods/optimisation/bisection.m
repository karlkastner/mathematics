% Wed  8 Aug 16:59:44 CEST 2018
% find threshold
%% bisection
function x0 = bisection(fun,xl,xr,n)
	fl = fun(xl);
	fr = fun(xr);
	if (fl==fr)
		error('not convex');
	end
	x0 = bisection_(fun,xl,xr,fl,fr,n)
end
function x0 = bisection_(fun,xl,xr,fl,fr,n)
	xc = 0.5*(xl+xr);
	if (n>0)
	fc = fun(xc);
	if (fc)
		% between fc and false side
		if (fl)
			x0 = bisection_(fun,xc,xr,fc,fr,n-1);
		else
			x0 = bisection_(fun,xl,xc,fl,fc,n-1);
		end
	else
		% between fc and true side
		if (fl)
			x0 = bisection_(fun,xl,xc,fl,fc,n-1);
		else
			x0 = bisection_(fun,xc,xr,fc,fr,n-1);
		end
	end
	else
		% interpolate final value to mid-section
		x0 = xc;
	end
end

