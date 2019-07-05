% Sun 12 Jun 13:43:27 CEST 2016
% Karl Kastner, Berlin
%
%% non-linear conjugate gradient
%% input:
%% x   : nx1 start vectort
%% opt : struct options
%% fdx : gradient constraint
function [x f r iter] = nlcg(afun,x,opt,lb,ub)
	if (nargin() < 3 || isempty(opt))
		opt = struct;
	end
	if (nargin() < 4)
		lb = -inf;
	end
	if (nargin() < 5)
		ub =  inf;
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
		opt.verbose = 0
	end
	if (~isfield(opt,'project'))
		opt.project = [];
	end
	if (isfield(opt,'ls_solver'))
		ls_solver = opt.ls_solver;
	else
		ls_solver = @line_search_quadratic2;
	end
	if (~isfield(opt,'ls_maxiter'))
		opt.ls_maxiter = 15;
	end
	if (~isfield(opt,'ls_h'))
		opt.ls_h = 1;
	end

	% calculate the residual at the initial step
%	[f r] = afun(x);
%	f  = sqrt(f'*f);
	f0 = [];
	p  = [];
	fold = inf;
	rold = zeros(size(x));
%	a  = 1;
%	p = r;
	iter = 1;
	while (1)
		% evaluate objective function and compute gradient at current point
		% r = r - a Ap
		% A'*residual = gradient : A'(Ax-b) = A'r
	
		% objective function value and Jacobien
		[f J] = afun(x);
		% gradient
		r = J'*f;
		% TODO distinguish between f and nf
		if (isvector(f))
			f = norm(f);
		end
		if (isempty(f0))
			% set initial search direction into direction of the gradient
			f0 = f;
			p  = r;
		else
		% update search direction
		%b = r'*r/(rold'*rold);
		dr = r-rold;
		% Polak-Riberie weight of last search direction
		b = r'*dr/(rold'*rold);
		% reset of search direction
		b = max(b,0);
		% ubdate search direction
		p = r + b*p;
		% reset search direction
		% This should not happen
		if (p'*r < 0)
			p = r;
			warning('reset of search direction');
		end

		end
		fold = f;
		
		% stop if the value of the objective function or gradient is invalid
		if (~isfinite(f) || any(~isfinite(r)))
			fprintf(['Stopped at iteration %d, without converging, ' ...
			       ,'as the objective function value or gradient was NaN\n'], ...
				 iter);
			break;
		end
	
		% perform line search in search direction and step
		% a = r'r / (p'Ap)
		% x = x + a*p;
		[x f a] = feval(ls_solver,afun,x,f,r,p,opt.ls_h,lb,ub,opt.ls_maxiter,opt.verbose,[]);

		if (~isfinite(f))
			break;
		end

		% status report
		if (opt.verbose)
		    fprintf('Iteration: %d Objective function value %g reduction %g\n', ...
									iter, f, fold-f);
			%fprintf('Objective function value: %f delta: %f \n',f,fold-f);
			fprintf('Gradient norm: L1 %g L2 %g Linf %g\n', ...
							norm(r,1),norm(r,2),norm(r,'inf'));
			fprintf('Gradient direction g''g_new/||g||||g_new||): %g\n',(r'*rold)/(norm(r)*norm(rold)));
		end
		if (opt.verbose > 1)
			fprintf('Location Gradient\n');
			fprintf('%8g %8g\n', [rvec(x);rvec(r)]);
		end
		% stop if there is no reduction of the objective function value
		if (f >= fold-opt.abstol)
			fprintf(['Stopped at iteration %d, without coverging, as there ' ...
				,'was no reduction of the objective function value\n'],iter);
			break;
		end
		% stop if objective function value sufficiently reduced (convergence possible)
		if (f < opt.reltol*f0)
			fprintf(['Stopped at iteration %d, as the value of the objective' ...
			       ,' function converged below the relative tolerance\n'],iter);
			break;
		end
		% stop if number of iterations exceeded
		if (iter >= opt.maxiter)
			fprintf(['Stopped at iteration %d, as maximum number of' ...
			       ,' iterations was reached.\n'],iter);
			break;
		end
		rold = r;
		iter = iter+1;
	end % while 1
end % nlcg

