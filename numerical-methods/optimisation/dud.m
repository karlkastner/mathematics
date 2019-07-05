% Mon 15 Aug 15:43:17 CEST 2016
% Karl Kastner, Berlin
%% optimization by the dud algorithm
% objective function must return a vector
% derivative free akgorithm for non-linear least squares
function [xnew nf] = dud(fun,x0,dx,opt,lb,ub)

	if (nargin() < 4 || isempty(opt))
		opt = struct;
	end
	if (nargin() < 5)
		lb = -inf;
		%fdx = true(size(x));
	end
	if (nargin() < 6)
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

	p = length(x0)
	f = NaN(p+1,1);
	X = NaN(p,p+1);
	X(1:p,1) = x0;
	X(end,:) = 1;
	if (1 == length(dx))
		dx = dx*ones(size(x0));
	end

	% vector valued objective function at initial point
	f0   = flat(fun(x0));
	nf0  = f0'*f0;
	nf   = nf0;
	n    = length(f0);
	f    = [f0,zeros(n,p)];

	% build up of the initial matrix
	for idx=2:p+1
		% solve problem in sub-space
		% TODO
		% perturb
		X(:,idx)     = x0;
		X(idx-1,idx) = x0(idx-1)+dx(idx-1);
		f(:,idx)     = flat(fun(X(:,idx)));
	end
	% TODO: sort ?
	iter = 0;
	dX = [];
	dF = [];
	while (1)
		nfold = nf;
		% determine the new point as minimum of the ls-mimisation
		% TODO shift and scale
		for idx=1:p
			dX(:,idx) = X(:,idx+1) - X(:,1);
			dF(:,idx) = f(:,idx+1) - f(:,1); 
		end
		y = 0;
		a  = (dF'*dF) \ (dF'*(y-f(:,1)));
		xnew = X(:,1) + dX*a;

		% determine coefficients (gradient estimate)
%		r = (X') \ f;
		% determine new zero
		%xnew = (f-a(end))./a;
%		A = X(1:p,1:p);
%		xnew = A \ (f-r(end));

		% check, that the new point is a minimizer
		% x_0 = a*dx + (1-a)*x_{1}, eq. 7 in Ralston 1978, c.f. armijo-wolfe
		maxiter = 50;
		[xnew fnew nf] = ls(fun,xnew,nfold,X(:,1),maxiter);

		% evaluate function at the new point
		%fnew = flat(fun(xnew));
		%nf = fnew'*fnew;
		% discard the oldest point and add new point
		% (TODO worst point?)
		X = [xnew, X(1:p,1:end-1)];
		f = [fnew, f(:,1:end-1)];
		% TODO, check, that the new set of points span the parameter space

		% status report
		if (opt.verbose)
			fprintf('Iteration: %d Objective function value %g reduction %g\n',iter,f,fold-f);
			%fprintf('Objective function value: %f delta: %f \n',f,fold-f);
			fprintf('Pseduo gradient norm: L1 %g L2 %g Linf %g\n',norm(r,1),norm(r,2),norm(r,'inf'));
		end
		% stop if the value of the objective function or gradient is invalid
		if (isnan(nf)) % || any(isnan(r)))
			fprintf('Stopped at iteration %d, without converging, as the objective function value or gradient was NaN\n',iter);
			break;
		end
		% stop if there is no reduction of the objective function value
		if (nf >= nfold-opt.abstol)
			fprintf('Stopped at iteration %d, without coverging, as there was no reduction of the objective function value\n',iter);
			break;
		end
		% stop if objective function value sufficiently reduced (convergence possible)
		if (nf < opt.reltol*nf0)
			fprintf('Stopped at iteration %d, as the value of the objective function converged below the relative tolerance\n',iter);
			break;
		end
		% stop if number of iterations exceeded
		if (iter >= opt.maxiter)
			fprintf('Stopped at iteration %d, as maximum number of iterations was reached.\n',iter);
			break;
		end
		iter = iter+1;
	end % while (1)
	% supply covariance estimate
	C = dX*inv(dF'*dF)*dX';
end % dud

function [x f nf] = ls(fun,x,nfold,xold,maxiter)
	iter = 0;
	x0 = x;
	while (true)
		%x_ = a*x + (1-a)*xold;
		f  = fun(x);
		nf = f'*f;
		[nf nfold]
		if (nf < nfold)
			iter
			break;
		end
		iter = iter+1;
		a = -((-0.5)^iter);
		x = a*x0 + (1-a)*xold;
		%x  = 0.5*(x + xold);
		%a=0.5*a;
		if (iter >= maxiter)
			printf('Maximum number of iteratios reached in ls\n');
			break;
		end
	end
end

