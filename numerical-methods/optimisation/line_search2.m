% Thu  9 Jun 20:02:35 CEST 2016
% Karl Kastner, Berlin
%
%% bisection method
%%
%% fun : objective funct
%% x0  : start value
%% f0  : objective function value at x0
%% g   : gradient at x0
%% p   : search direction from x0 (p = g for steepest descend)
%% h   : initial step length (default 1)
%% lb  : lower bound for x
%% up  : upper bound for x
function [x f h cflag] = line_search2(fun,x0,f0,g,p,h,lb,ub,maxiter);
	if (isempty(p))
		p = g;
	end
	if (isempty(h))
		h = 1;
	end
	if (nargin() < 7)
		lb = -inf;
	end
	if (nargin() < 8)
		ub = inf;
	end
	if (nargin() < 9 || isempty(maxiter))
		maxiter = 12;
	end
	% ramp up step length
	% h = 2*h;
	cflag = false;

	% evaluate function at left
	hl = 0;
	fl = f0;

	% function value at right end of interval
	hr = h;
	xr = max(min(x0-hr*g,ub),lb);
	fr = fun(xr);
	fr = norm(fr);

	iter = 0;
	while (1)
		% check armijo-wolfe conditions for sufficient step length
		if (1)
		if (fl < fr)
			if (armijo_stopping_criterion(f0,g,fl,hl))
				f = fl;
				h = hl;
				cflag = true;
				break;
			end
		else
			if (armijo_stopping_criterion(f0,g,fr,hr))
				f = fr;
				h = hr;
				cflag = true;
				break;
			end
		end
		end % if (1)

		% update central value
		hc = 0.5*(hl+hr);
		xc = max(min(x0-hc*g,ub),lb);
		fc = fun(x0-hc*g);
		fc = norm(fc);
	
		% discard larger half
		if (fl > fr)
			% discard left half
			fl = fc;
			hl = hc;
		else
			% discard right half
			fr = fc;
			hr = hc;
		end
		iter = iter+1;
		if (iter >= maxiter)
			if (fl < fr)
				f = fl;
				h = hl;
			else
				f = fr;
				h = hr;
				
			end
			break;
		end
	end %  while (1)
	x = max(min(x0 - h*g,ub),lb);
	% x_ = [0; 0.5; 1]; poly = PolyOLS(2); p=poly.regress(x_,[9.637482791709112e+03; 8.346538576065512e+03; 7.890249687264490e+03]); dp = [p(2);2*p(3)], xz=-dp(1)/dp(2)
end % line_search2

