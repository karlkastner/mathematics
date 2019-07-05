% Sun 21 Aug 16:27:27 CEST 2016
% Karl Kastner, Berlin
%
%% quadratic line search
%% fun : objective funct
%% x0  : start value
%% f0  : objective function value at x0
%% g   : gradient at x0
%% dir : search direction from x0 (p = g for steepest descend)
%% h   : initial step length (default 1)
%% lb  : lower bound for x
%% up  : upper bound for x
function [x f dx] = line_search_quadratic(fun,x,f0,g,dir,h,lb,ub,maxiter,verbose)
	verbose = 1;
	if(isempty(h))
		h = 1;
	end
	iter = 0;
	while (true)
		% get projected hessian and gradient
		% TODO check if the function returns hessian and gradient
		[Hp gp f] = hessian_projected(fun,x,dir,h);
		fprintf('Line search iteration: %d dx: %g f: %g\n',iter,0,f);
		% optimal step length along gradient
		dx = Hp \ gp;
		while (true)
			% step
			xnew = x - dx*dir;
			% apply constraints
			xnew = max(min(xnew,ub),lb);
			% evaluate function at new point
			fnew = fun(xnew);
			fprintf('Line search backtrace dx: %g f: %g\n',dx,fnew);
			if (~isfinite(fnew))
				%error('line search');
				break;
			end
			% TODO check armillo-wolfe
			if (fnew < f0)
				break;
			end
			% reduce step length
			dx = 0.5*dx;
			% TODO limit number of iterations
		end
		x = xnew;
		f = fnew;
		if (~isfinite(f))
			x = NaN*x;
			break;
		end
		% check for convergence
		iter = iter+1;
		if (iter >= maxiter)
			break;
		end
	end
end % line_search_quadratic

