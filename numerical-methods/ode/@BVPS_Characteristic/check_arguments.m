% Sat 28 Oct 14:43:06 CEST 2017
function bvp2check_arguments(obj)
%	if (nargin()<1)
%		obj.opt = struct();
%	end

	% number of grid points
	if (~isfield(obj.opt,'nx'))
		obj.opt.nx = 100;
	end

	if (~isfield(obj.opt,'xs'))
		obj.opt.xs = 1;
	end

	% obj.options for solver
	if (~isfield(obj.opt,'sopt'))
		obj.opt.sopt = struct();
	end

	% maximum number of iterations
	if (~isfield(obj.opt.sopt,'maxiter'))
		obj.opt.sopt.maxiter = opt.nx;
	end
	% relaxation constant
	if (~isfield(obj.opt.sopt,'relaxation'))
		obj.opt.sopt.relaxation = 0.5;
	end
	
	if (~isfield(obj.opt,'reconstruct_y'))
		obj.opt.reconstruct_y = true;
	end

	if (~isfield(obj.opt,'bcarg'))
		obj.opt.bcarg = {};
	end

	if (~isfield(obj.opt,'balance'))
		obj.opt.balance = false;
	end

	% minimum number of grid points
	obj.opt.nx = max(2,obj.opt.nx);
end 
