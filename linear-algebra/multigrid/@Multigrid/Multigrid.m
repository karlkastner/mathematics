% 2024-01-10 19:59:21.792153268 +0100
% Karl Kastner, Berlin
%
%% multigrid on structured 2D grid
classdef Multigrid < handle
	properties
		% a : nvar x nvar x 1 or nvar x nvar x (nx*ny)
		% coefficients of the reaction-part 
		a
		% ad : 2 x 3 x nvar : coefficient of the advection part
		ad
		% e : 2 x nvar : coefficiets of the diffusion part
		e
		% advection-diffusion kernel
		% kernel = struct('ad',[]);
		% extend of domain
		L
		% number of grid cells
		n
		% jacobi-iteration parameters
		o = 2/3;
		% level coefficients and variables
		s = struct();
		% number of independent variables
		nvar
		% spatial discretization dx = L/n
		dx
		fun
		% 
		% opt = struct('ndiag',1);
		opt = struct('operatord_dependent_grid_transfer', false, ...
			     'downsampling_mode','fw' ...	
				... % relative tolerance for stopping iterations
				,'reltol', sqrt(eps) ...
				,'abstol', 1e-12 ...
				... % maximum number of cycles
				,'maxiter', 1000 ...
				... % number of pre and post-smoothing iterations per cycle
				,'npresmooth',  1 ...
				,'npostsmooth', 1 ...
				... % cycle counter, 1 = v-cycle, 2 = w-cycle
				,'gamma', 1 ...
				);
	end % properties
	methods
		function obj = Multigrid()
		end % Multigrid
	end % method
end % classdef Multigrid

