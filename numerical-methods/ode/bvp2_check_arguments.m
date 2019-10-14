% Sat 28 Oct 14:43:06 CEST 2017
function opt = bvp2check_arguments(opt)
	if (nargin()<1)
		opt = struct();
	end

	% number of grid points
	if (~isfield(opt,'nx'))
		opt.nx = 100;
	end

	if (~isfield(opt,'xs'))
		opt.xs = 1;
	end

	% options for solver
	if (~isfield(opt,'sopt'))
		opt.sopt = struct();
	end

	% maximum number of iterations
	if (~isfield(opt.sopt,'maxiter'))
		opt.sopt.maxiter = opt.nx;
	end
	% relaxation constant
	if (~isfield(opt.sopt,'relaxation'))
		opt.sopt.relaxation = 0.5;
	end

	% minimum number of grid points
	opt.nx = max(2,opt.nx);
end 
