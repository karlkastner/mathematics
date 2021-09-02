% Tue 15 Jun 19:42:20 CEST 2021
% c.f. walker & Li 2011 
% note that this does not perform better than picard iteration
function [x_k, cflag, iter] = anderson_iteration(fun,x0,opt)
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
	if (nargin()<4)
		flag = true;
	end
	iter = 0;
	cflag = 0;
	x_k = fun(x0);
	X = [];
	F = [];
	while (1)
		iter = iter+1;
		X   = [X,x_k-x0];
		ndx = norm(X(:,end));
		% check for convergence
		if (ndx < opt.reltol*norm(x_k))
			cflag = 1;
			break;
		end

		%if (iter > 1 && ndx < opt.reltol*ndxold)
		%	cflag = 2;
		%	break;
		%end
		if (iter > opt.maxiter)
			warning('stopped before converging because maximum number of iterations was reached');
			break;
		end

		f_k   = fun(x_k) - x_k;
		F     = [F,f_k];
		gamma = F \ f_k;
		x0    = x_k;
		%x_k   = x_k + f_k - (X*f_k - F*f_k)*gamma;
		x_k   = x_k + f_k - (X + F)*gamma;
	end % while
end % anderson_iteration

