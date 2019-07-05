% Mon 22 Aug 12:03:03 CEST 2016
% Karl Kastner, Berlin
%
%% polynomial line search
%% fun : objective funct
%% x0  : start value
%% f0  : objective function value at x0
%% g   : gradient at x0
%% dir : search direction from x0 (p = g for steepest descend)
%% h   : initial step length (default 1)
%% lb  : lower bound for x
%% up  : upper bound for x
function [x f dx] = line_search_polynomial(fun,x,f0,g,dir,h,lb,ub,maxiter,order,verbose)
	verbose = 1;
	if (nargin() < 10 || isempty(order))
		order = 4;
	end
	dir = dir/norm(dir);
	btmax = 5;
	iter  = 0;
	while (true)
		% get projected hessian and gradient
		% TODO check if the function returns numerical hessian and gradient
        	% get fourth projected derivative and series expansion
		[d p] = directional_derivative(fun,x,dir,h,order);
		f = p(end);
		% get extrema of series expanded function
		[dx y minflag realflag] = poly_extrema(p);
        	% TODO check if minima are real                                             

		% select minimum
		dx_ = dx(minflag);
		switch (length(dx_))
		case {0}
			error('here');
		case {1} % 
			dx = dx_;
		otherwise % choose the minima closest to the old point
			dis = abs(dx_);
			[mdis mdx] = min(dis);
			dx = dx_(mdx);
		end
		% invert to obtain optimal step length along gradient
		dx = -dx;

		fprintf('Line search iteration: %d dx: %g f: %g\n',iter,0,f);
		btx = 0;
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
			% TODO limit number of iterations
			btx = btx+1
			if (btx > btmax)
				break;
			end
			% reduce step length
			dx = 0.5*dx;
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
end % line_search_quartic

