% 2016-05-08 10:56:59.854888599 +0200
%% rectangular window
function win = rectwin(x,x0,L)
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
		win = abs(x-x0) < 0.5*L;
		win = win/sum(win);
	end	
end

