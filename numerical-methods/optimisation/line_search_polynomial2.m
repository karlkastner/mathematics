% Thu  1 Sep 18:09:08 CEST 2016
% Karl Kastner, Berlin
%
%% cubic line search
%% fun : objective funct
%% x0  : start value
%% f0  : objective function value at x0
%% g   : gradient at x0
%% dir : search direction from x0 (p = g for steepest descend)
%% h   : initial step length (default 1)
%% lb  : lower bound for x
%% up  : upper bound for x
function [x f dx] = line_search_cubic2(fun,x,f0,g,dir,h,lb,ub,maxiter,order,verbose)
	verbose = 1;
	dir   = dir/norm(dir);
	btmax = 5;
	iter  = 0;
	
	a = 0;
	f = f0;
	
	% initia linear step
	a(2) = g'*d/(d'*d);
	f(2) = fun(a(2));

	% initial quadratic step
	p    = polyfitd(a,f,a(1),g'*d/sqrt(d'*d));
	p1   = polyder(p);
	a(3) = roots(p1);
	f(3) = fun(a(3));

	% initial cubic step
	p = polyfitd(a,f,a(1),g'*d/sqrt(d'*d));
	% case differentiation here necessary
	p = 


	if (2*an < a(2)+a(3))
		% drop a4
		if (an < a(1))
			a = [an,a(1:3)];
		else if (an < a(2))
			a = [a(1), an, a(2:3)];
		else
			a = [a(1:2),an,a(3)];
		end
	else
		% drop a1
		if (an > a(4))
			a = [a(2:4),an];
		else if (an > a(3))
			a = [a(2:3), an, a(4)];
		else
			a = [a(2),an,a(3:4)];
		end
	end

	
	% iteration
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

