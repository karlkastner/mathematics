% 2016-05-08 10:56:59.854888599 +0200
%% rectangular window
function [win,Lw] = rectwin(x,x0,L)
	if (isscalar(x))
		n = x;
		win = ones(n,1)/n;
	else
		if (nargin() < 2)
			x0 = 0.5*(x(1)+x(end));
		end
		if (nargin() < 3)
			L = x(end)-x(1);
			% scale
			n = length(x);
			L = L*n/(n-1);
		end
		n     = length(x);
		dx    = (x(end)-x(1))/(n-1);
		win1  = abs(x-x0) < 0.5*L;
		L1    = dx*sum(win1);
		if (L1>L)
			win2  = abs(x-x0) < 0.5*L-dx;
		else
			win2  = abs(x-x0) < 0.5*L+dx;
		end
		L2    = dx*sum(win2);
		p     = (L-L2)/(L1-L2);
		win   = p*win1 + (1-p)*win2;
		Lw    = sum(win)*dx;
		win = win/sum(win);
	end	
end

