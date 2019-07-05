% 2016-08-28 23:16:58.134778707 +0200
% Karl Kastner, Berlin
%% line search by wolfe method
%% c.f.: OPTIMIZATION THEORY AND METHODS - Nonlinear Programming, Sun, Yuan
function [x f a cflag iter] = line_sear_wolfe(fun,x0,f0,g,p,a0,lb,ub,maxiter,verbose)
	if (nargin() < 10)
		verbose = false;
	end
	% initialise
	al = 0;
	a  = a0;
	a = g'*p/(p'*p);
	ar = inf;
	cflag = false;
	t    = 2;
	rho  = 0.0;% must be < 1/2 for superlinear convergence
	iter = 0;
	xold=x0;
	fold=f0;
	aold=a0;

	while (1)
		x = x0-a*p;
		f = fun(x);
		if (verbose)
			fprintf('LS iter: %d a: %g f: %g liml: %g limr %g\n',iter,a,f0-rho*a*g'*p,f,f0-(1-rho)*a*g'*p);
		end
%		iter
%		disp(x0-[al a ar]*p)
%		disp([f0-rho*a*g'*p f f0-(1-rho)*a*g*p])
%		disp(f<[f0-rho*a*g'*p f0 f0-(1-rho)*a*g*p])
%		pause

		% note, this is smaller equal in tsu, but has to be strictly smaller
		% note the wrong sign in tsu in both equations
		% 2.5.3
		if (f >= f0 - rho*a*g'*p)
			% no reduction at this distance
			% step size too large
			% reduce step size
			ar=a;
			a = 0.5*(al+ar);
		else
			% 2.5.4
			if (f < f0 - (1-rho)*a*g'*p )
				% reduction too low
				% step size too small
				% increase step size
				al = a;
				if (isfinite(ar))
					a = 0.5*(al+ar);
				else
					a = t*a;
				end
			else
				% reduction sufficient, stop
				cflag = true;
				% this can happen if rhs = inf
				if (fold < f)
					a=aold;
					x=xold;
					f=fold;	
				end
				break;
	% if rho ~= 0 : check if boundaries are better:
	%			if (fl < f)
	%				a = al;
	%				f = fl;
	%			end
	%			if (fr < f)
	%				a = ar;
	%				f = fr;
	%			end
			end
		end 
		iter = iter+1;
		if (iter > maxiter)
			break;
		end

	xold=x;
	fold=f;
	aold=a;
	end % while
end % line_search_wolfe

