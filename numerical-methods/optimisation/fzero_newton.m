% Fri 24 Jan 09:42:23 +08 2020
function x = fzero_newton(fun,grad_f,x,opt,LB)
	if (nargin()<4)
		opt = optimset();
	end
	if (isempty(opt.MaxIter))
		% to be sqrt(eps) precise
		opt.MaxIter = round(log2(1/sqrt(eps)));
	end
	k = 0;
	while (1)
		k = k+1;
		f = fun(x);
		g = grad_f(x);
		x = x - (f./g);
		if (nargin()>4)
		x = max(LB,x);
		end

		if (k > opt.MaxIter)
			break;
		end
	end % while 
end

