% Wed  8 Jun 12:05:58 CEST 2016
% Karl Kastner, Berlin
%
%% solve non-linear system by cuachy's method
%% slower than quadratic optimisation, but does not require a hessian
%
%% fun : objective function, returns
%%	f  : scalar, objective function value 
%%	g  : nx1, gradient
%% x   : nx1, initial position
%% opt : options
function [x f g] = cauchy2(fun,x,opt,fdx)
	if (nargin() < 3)
		opt = struct;
	end
	if (~isfield(opt,'maxiter'))
		opt.maxiter = 100;
	end
	if (~isfield(opt,'abstol'))
		opt.abstol = sqrt(eps);
	end
	if (~isfield(opt,'reltol'))
		opt.reltol = 1e-3;
	end
	if (~isfield(opt,'verbose'))
		opt.verbose = false;
	end
	if (~isfield(opt,'project'))
		opt.project = [];
	end
	if (isfield(opt,'ls_solver'))
		ls_solver = opt.ls_solver;
	else
		% line_search2
		ls_solver = @line_search_quadratic2;
	end
	if (~isfield(opt,'ls_maxiter'))
		opt.ls_maxiter = 15;
	end
	fold  = inf;
	[f g] = fun(x);
	% initial step length
	h = 1./sqrt(g'*g);

	f = sum(f);
	f0     = f;
	iter   = 0;
	while (1)
		if (opt.verbose)
			fprintf('Iteration: %d Objective function value %1.15f\n',iter,f);
		end
		% stop if the value of the objective function is invalid
		if (isnan(f))
			fprintf('Stopped at iteration %d, without converging, as the objective function value was NaN\n',iter);
			break;
		end
		% stop if there is no reduction of the objective function value
		if (f >= fold-opt.abstol)
			fprintf('Stopped at iteration %d, without coverging, as there was no reduction of the objective function value\n',iter);
			break;
		end
		if (f < opt.reltol*f0)
			fprintf('Stopped at iteration %d, as the value of the objective function converged below the relative tolerance\n',iter);
			break;
		end

		% save function value
		fold = f;

		% step
		xold    = x;
		if (nargin() > 3)
			g(~fdx) = 0;
		end
		%/home/pia/phd/src/lib/optimisation/line_search2.m
                [x f h] = ls_solver(fun,x,f,g,g,h,-inf,inf,opt.ls_maxiter,1,opt);
		f = sum(f);
		
%		if (~isempty(opt.project))
%			x = feval(opt.project,x,xold);
%			% recompute function value
%			f = fun(x);
%			f = sum(f);
%		end
		
		% increase interation counter
		iter = iter+1;
		if (iter >= opt.maxiter)
			fprintf('Stopped at iteration %d, as maximum number of iterations was reached.\n',iter);
			break;
		end
		% evaluate objective function at new location
		[f g] = fun(x);
		f = sum(f);
	end % while
end % function cauchy2

