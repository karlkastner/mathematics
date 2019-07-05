% Mon  5 Sep 00:15:50 CEST 2016
% Karl Kastner, Berlin
%% non-linear least squares
function [x res nres g ng A] = nlls(fun,x,opt,lb,ub)
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
	if (~isfield(opt,'logfile'))
		opt.logfile = [];
	end
%	if (~isfield(opt,'project'))
%		opt.project = [];
%	end
	if (isfield(opt,'ls_solver'))
		ls_solver = opt.ls_solver;
	else
		ls_solver = @line_search_quadratic2;
	end
	if (~isfield(opt,'ls_maxiter'))
		opt.ls_maxiter = 15;
	end
	if (~isfield(opt,'lambda'))
		opt.lambda = 1e-1;
	end
	% TODO, no magic numbers
        %reltol = 1e-3;
	bt = false;

	nlls = struct();
	nlls.x       = x;

%	res0 = [];
	% residual norms [intial, 1, 2, 3 ...]
	nres = [];
%	an   = [];
	% gradient of last step
	gold = NaN(size(x));
	iter = 0;
	while (1)
		iter = iter+1;
		xold = x;

		% residual and Jacobien
		[res A] = fun(x);
		if (1 == iter)
			nres(1) = norm(res);
			%nresold;
			%nresold = norm(res);
			%nres0 = nres(1);
		end
		%[nresold nres(iter)]

		% update function value
		% TODO le-ma stabilistation if A singular
		%dx = (A \ res);
		% make A'A diagonally dominant,
		% this guarantees that the search direction is a descend direction
		AA = A'*A;
	        DD = diag(diag(AA));
	        %[lmin lmax ratio] = gershgorin_shift(A'*A);
		%ratio = max(ratio);
		%if (ratio < 0)
		%	ratio = 0;
		%end
		% lambda = ratio
		%D = diag(diag(AA);
		% TODO, use RRQR
		iter_ = 0;
		% gradient
		g      = A'*res;
		% the singular values of AA
		% (TODO can be used for inversion)
		s      = eig(AA);
		maxs   = max(s);
		I      = eye(size(AA));

		% TODO, update gradient only in the direction of the largest sv
		% maybe it is smarter not to update maxs during the iteration
		while (1)
			%dx = (AA + opt.lambda*DD)\g;
			dx           = (AA + opt.lambda*maxs*I) \ g;
			% step
			x            = xold - dx;
			res          = fun(x);
			nres(iter+1) = norm(res);
			if (nres(iter+1) < nres(iter))
				opt.lambda = 0.5*opt.lambda;
				break;
			else
				fprintf(1,'No reduction in function value, increasing lambda\n');
				opt.lambda  = 2*opt.lambda;
			end
			iter_ = iter_+1;
			if (iter_ > opt.ls_maxiter)
				fprintf(1,'Iteration stopped, as there was no reduction in function value\n')
				ng(iter) = NaN;
				return
			end
		end % while (1)

		% gradient norm
		ng(iter,1) = norm(g);

		% store initial residual
		%if (isempty(res0))
		%	ng0   = ng(1);
		%	nlls.nres(1) = nres(1)
		%end

		% status report
		if (opt.verbose)
		    fprintf('Iteration: %d Objective function value %g reduction %g\n', ...
									iter, nres(end), nres(end-1)-nres(end));
			%fprintf('Objective function value: %f delta: %f \n',f,fold-f);
			fprintf('Gradient norm: L1 %g L2 %g Linf %g\n', ...
							norm(g,1),norm(g,2),norm(g,'inf'));
			fprintf('Change of gradient: g''g_old/||g||||g_old||): %g\n',(g'*gold)/(norm(g)*norm(gold)));
			fprintf('Difference between search direction and gradient: p''g/||p||||g||): %g\n',(dx'*g)/(norm(dx)*norm(g)));
		end
		if (opt.verbose > 1)
			fprintf('Location Gradient\n');
			fprintf('%8g %8g\n', [rvec(x);rvec(g)]);
		end
		% log
		if (length(opt.logfile) > 0)
			% gradient
			nlls.g(:,iter)    = g;
			% step (search direction times step length)
			nlls.dx(:,iter)   = dx;
			nlls.x(:,iter+1)  = x;
			nlls.lambda(iter) = opt.lambda;
			nlls.ng(iter)     = ng(iter);
			nlls.nres(iter+1) = nres(iter+1);
			nlls.A(:,:,iter)  = A;
			save(opt.logfile,'nlls');
		end

%		% stop if there is no reduction of the objective function value
%		if (f >= fold-opt.abstol)
%			fprintf(['Stopped at iteration %d, without coverging, as there ' ...
%				,'was no reduction of the objective function value\n'],iter);
%			break;
%		end
		% stop if objective function value was sufficiently reduced (convergence possible)
		if (nres(end) < opt.reltol*nres(1)) %nres0)
			fprintf(['Stopped at iteration %d, as the value of the objective' ...
			       ,' function converged below the relative tolerance\n'],iter);
			break;
		end
		% stop if gradient is sufficiently small (convergence)
		if (ng(end) < opt.reltol*norm(x) || ng(end) < opt.abstol*sqrt(length(g)))
			fprintf(1,['Stopped at iteration %d, as the gradient of the objective' ...
				  ,' function converged below the relative tolerance\n']);
			break;
		end
		% stop if number of iterations exceeded
		if (iter >= opt.maxiter)
			fprintf(['Stopped at iteration %d, as maximum number of' ...
			       ,' iterations was reached.\n'],iter);
			break;
		end
%		rold = r;
%		iter = iter+1;

		gold    = g;
		%nresold = nres;
	end % while 1
end % nlls

