% Tue 16 Aug 17:22:59 CEST 2016
% Karl Kastner, Berlin
%
%% least squares by the secant method
%% Barnes, 1965
%% Wolfe, 1959
%% Fletcher 1980, 6.3
%% seber 2003
%% gerber

function [x resn r J] = ls_generalized_secant(afun,x,opt,lb,up)
	if (nargin() < 3 || isempty(opt))
		opt = struct;
	end
	if (nargin() < 4)
		lb = -inf;
		%fdx = true(size(x));
	end
	if (nargin() < 5)
		ub =  inf;
		%fdx = true(size(x));
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
	if (~isfield(opt,'lsmaxiter'))
		opt.lsmaxiter = 0;
	end
	% initial value
	r = afun(x);
	f0 = norm(r);
	n = length(r);
	p = length(x);
	% initial jacobian
	J = rand(length(r),p);
'svd'
	svd(J)
	%[length(eye(p);zeros(
	D = zeros(p,p);
	iter = 0;
	f = NaN;
	r = 0;
	while (1)
		xold = x;
		fold = f;
		rold = r;

		% get the residual (r = Ax-b), but not the gradient
		r    = afun(x);
		resn = norm(r);
		f    = resn;
		% step
		%xnew = x + J \ r;
		xnew = x - (J'*J) \ J'*r; % +?
		[x f h cflag] = line_search2(afun,xnew,f,xnew-x,xnew-x,[],-inf,inf,12);
		% TODO line search
		x = xnew;
		y  = r - rold;
		dx = x - xold;
		% the column space is only of length p, hence the oldest vector is always forgotten
		idx = mod(iter,p)+1;
		D(:,idx) = dx;
		% determine a vector that is orthogonal to previous steps
		w = ones(p,1);
		if (iter > 1)
			D_ = D;
			D_(:,idx) = 0;
			w = w - D*(D'*w);
		end
		w = w/norm(w);

		% update the Jacobian
		% TODO use udate of the inverse
		% NOTE had to be transposed to match size
J_ = J;
		J  = J + (1/(w'*dx))*(w*(y - J*dx)')';
norm(J-J_)
[x xold]
pause

		% status report
		if (opt.verbose)
		    fprintf('Iteration: %d Objective function value %g reduction %g\n' ...
									, iter, f, fold-f);
			%fprintf('Objective function value: %f delta: %f \n',f,fold-f);
			fprintf('Gradient norm: L1 %g L2 %g Linf %g\n' ...
							, norm(r,1),norm(r,2),norm(r,'inf'));
		end
		% stop if the value of the objective function or gradient is invalid
		if (isnan(f) || any(isnan(r)))
			fprintf(['Stopped at iteration %d, without converging, ' ...
			       ,'as the objective function value or gradient was NaN\n'] ...
				, iter);
			break;
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
		iter = iter+1;
	end % while 1
end % ls_generalized_secant


