% Sat 28 Oct 14:43:06 CEST 2017
function arguments(obj)

	% number of grid points
	if (isempty(obj.nx))
		obj.nx = 100;
	else
		% minimum number of grid points
		obj.nx = max(2,obj.nx);
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
		obj.opt.sopt.maxiter = sum(obj.nx);
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

	if (~isfield(obj.opt,'couple_frequency_components'))
		obj.opt.couple_frequency_components = true;
	end

end % check_arguments

