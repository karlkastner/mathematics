% Fri 24 Jan 19:24:37 +08 2020
% fails on x.^(2n) roots!
function x0 = fzero_bisect(fun,xl,xr,opt)
	if (nargin()<4 || ~isfield(opt,'MaxIter') || isempty(opt.MaxIter))
		MaxIter = 15;
	else
		MaxIter = opt.MaxIter;
	end
	yl = feval(fun,xl);
	yr = feval(fun,xr);
	x0 = fzero_bisect_(fun,xl,xr,yl,yr,MaxIter);
end

function xc = fzero_bisect_(fun,xl,xr,yl,yr,k)
		xc  = 0.5*(xl+xr);
		yc  = fun(xc);
%		[xl,xc,xr,yl,yc,yr]
		fdx = sign(yl) == sign(yc);
	
		xl(fdx) = xc(fdx);
		yl(fdx) = yc(fdx);
		%fdx = sign(yr) == sign(yc);
		fdx     = ~fdx;
		xr(fdx) = xc(fdx);
		yr(fdx) = yc(fdx);
%		[xl,xc,xr,yl,yc,yr]
%pause
	if (k>0)
		xc = fzero_bisect_(fun,xl,xr,yl,yr,k-1);
	end
end

