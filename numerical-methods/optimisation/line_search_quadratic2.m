% Wed 31 Aug 14:52:46 CEST 2016
%
%% quadratic line search
% This can be done smarter, at initial step, gradient is known, exploit this to step just with 2 points + 1 gradient
function [xn fn an cflag iter] = line_search_quadratic_2(fun,x0,f0,g,d,h,lb,ub,maxiter,verbose,opt)
	if (nargin() < 11)
		opt = struct();
	end
	

	tol = 0.01;
	h = 1;

	while (1)
		a = h*[0 1 2]*(d'*g)/(d'*d);
		f = [sum(f0),sum(fun(x0-a(2)*d)),sum(fun(x0-a(3)*d))];
		if (all(isfinite(f)))
			break;
		end
		h = 0.5*h;
	end
	%f = [f0,fun(x0-a(2)*d),fun(x0-a(3)*d)];
	
	aold = a(2);
	xold = x0;
	fold = f(2);
	
	% initial gradient for stopping
	% for the stopping criteria the gradient of the polynomial is used,
	% rather than the projected gradient
	p   = polyfit(a,f,2);
	p1  = polyder(p);
	gp0 = polyval(p1,0);

	
	cflag = false;
	iter = 0;
	while (1)

	% get new optimum
	an  = roots(p1);
	if (isempty(an))
		xn = xold;
		fn = fold;
		an = aold;
		clfag = false;
		warning('no minimum');
		break;
	end
	
	% evaluate function
	xn = x0 - an*d;
	if (isfield(opt,'project') && ~isempty(opt.project) && opt.project)
		xn = opt.project(x0,xn);
	end
	while (1)
		fn = sum(fun(xn));
		if (isfinite(fn))
			break;
		end
		xn = 0.5*(xn+xold);
	end

if (0)
	% backtrack to ensure reduction
	while (fn > f(2))
		%an = 0.5*(an+aold);
		an = 0.5*(an+a(2));
		xn  = x0 - an*d;
		fn  = sum(fun(xn));
%		% TODO limit iterations
	end
end
	
	% interleave new point
	if (an > a(2))
%		'right'
	if (0)
		if (an > a(3))
			a = [a(2:3), an];
			f = [f(2:3), fn];
		else
			a = [a(2), an, a(3)];
			f = [f(2), fn, f(3)];
		end
	else
		if (fn > f(2))
			if (fn > f(3))
				warning('no decrease in function value');
%				error('here')
			end
			a = [a(1:2) an];
			f = [f(1:2) fn];
		else
			if (fn > f(1))
				warning('no decrease in function value');
%				error('here')
			end
			a = [a(2) an a(3)];
			f = [f(2) fn f(3)];
		end
	end
	else
%	'left'
	if (0)
		if (an < a(1))
			a = [an, a(1:2)];
			f = [fn, f(1:2)];
		else
			a = [a(1), an, a(2)];
			f = [f(1), fn, f(2)];
		end
	else

		if (fn > f(2))
			if (fn > f(1))
				warning('no decrease in function value');
%				error('here')
			end
			a  = [an a(2:3)];
			f  = [fn f(2:3)];
		else
			if (fn > f(3))
				warning('no decrease in function value');
%				error('here')
			end
			a = [a(1) an a(2)];
			f = [f(1) fn f(2)];
		end
	end
	end

	% check for singularity
	if (    abs((a(2)-a(1))) < (eps).^0.25*abs(a(3)-a(2)) ...
	     || abs((a(3)-a(2))) < (eps).^0.25*abs(a(2)-a(1)) )
		disp('ls stopping as dx too small');
		break;
	end

	% determine cubic gradient, stop if sufficiently small
	% TODO, cubic may be insufficient for even functions
	%p  = polyfit([x x4]);
	%g  = polyval(p,xn);
	
	p  = polyfit(a,f,2);
	p1 = polyder(p);
	gp = polyval(p1,aold);

	%p  = polyfit();
	%p1 = polyder(p);
	% terminate if the line search sufficiently exact
	%if (abs(fn-fn_) <= tol*(fn-f0))
	%end
	%p1 = polyfit([xn xold],[fn fold]);
	%p0 = polyder(p);
	%g  = p0; % == polyval(p0,x), as constant
	if (verbose)
		printf('Line search iteration %d gradient previous iteration, f: % g: %g g/g0: %g\n',iter,fn,gp,gp/gp0);
	end

%	[abs(gp),abs(gp0)]
%pause	
	if (abs(gp) < tol*abs(gp0) && fn < f0)
		cflag = true;
		break;
	end

	iter = iter + 1;
	if (iter > maxiter)
		disp('ls stopping as maxiter reached');
		break;
	end
	
	aold  = an;
	fold = fn;
	
	end % while

end % function

