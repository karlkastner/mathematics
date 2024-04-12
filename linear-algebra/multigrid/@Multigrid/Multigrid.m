% 2024-01-10 19:59:21.792153268 +0100
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
		% cycle counter, 1 = v-cycle, 2 = w-cycle
		gamma = 1;
		% number of pre and post-smoothing iterations per kevek
		m = 1;
		% maximum number of cycles
		maxiter = 1000;
		% relative tolerance for stopping iterations
		reltol = sqrt(eps);
		% level coefficients and variables
		s = struct();
		% number of independent variables
		nvar
		% spatial discretization dx = L/n
		dx
		% 
		% opt = struct('ndiag',1);	
	end % properties
	methods
		function obj = Multigrid()
		end % Multigrid
	end % method
end % classdef Multigrid

