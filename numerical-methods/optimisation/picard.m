% Fri  8 Sep 08:59:04 CEST 2017
%% picard iteration
function [x cflag] = picard(fun, x, opt)
	if (nargin()<3)
		opt = struct();
	end
	if (~isfield(opt,'abstol'))
		opt.abstol  = 1e-7;
	end
	if (~isfield(opt,'reltol'))
		opt.reltol = 1e-3;
	end
	if (~isfield(opt,'maxiter'))
		opt.maxiter = length(x);
	end
	if (~isfield(opt,'relaxation'))
		opt.relaxation = 1;
	end
	if (~isfield(opt,'verbose'))
		opt.verbose = 0;
	end
	if (~isfield(opt,'miniter'))
		opt.miniter = 0;
	end

	ndx = inf;
	iter = 0;
	cflag = false;
	while (1)
		iter   = iter+1;
		ndxold = ndx;
		xold   = x;

		% new function value
		x      = fun(x);

		% step direction
		dx     = x-xold;

		% step length
		ndx   = norm(dx);

		% relaxation
		% with newtons method, the number of significant digits changes
		% p = min(1, ndx/ndxold);

		% update (relaxed step)
		x     = xold + opt.relaxation*dx;
		% dxold = 

%		clf
%		subplot(2,2,1)
%		plot(select(real([xold]),1:100))
%		hold on
%		plot(select(real([x]),1:100))
%		subplot(2,2,2)
%		x_ = select(x,101:200);
%		plot([imag(x_),real(x_)])
%		pause

		switch (opt.verbose)
		case {1}
			fprintf('||x-xold|| %g\n',ndx);
		case {2}
			fprintf('||x-xold|| %g\n',ndx);
			disp(x);
		end

		if (iter >= opt.miniter)		
		if (ndx < opt.reltol*norm(x))
			cflag = 1;
			break;
		end

		if (iter > 1 && ndx < opt.reltol*ndxold)
			cflag = 2;
			break;
		end
		if (iter > opt.maxiter)
			warning('stopped before converging because maximum number of iterations was reached');
			break;
		end
		end
	end % while
end % function picard

