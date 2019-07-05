% Thu Apr 28 04:12:31 MSD 2011
% Sat Feb 19 19:59:39 MSK 2011
% 2013/04/14 23:45
% Sa 6. Feb 23:17:44 CET 2016
% Karl Kästner, KTH Stockholm, Berlin
%
%% Lax-Friedrich-Method
%% for hyperbolic conservation laws
%% err = O(dt) + O(dx)
%% |a dt/dx| < 1
classdef Lax_Friedrich < Finite_Volume
	properties
		% averaging matrices
		M
		% difference matrices
		Dc
	end % properties
	methods
		function obj = init(obj, X, n)
			obj = init@Finite_Volume(obj, X,n);

			% central average matrix
			M = spdiags(ones(n,1)*[1 0 1], -1:1, n, n);
			% central difference matrix
			D = spdiags(ones(n,1)*[-1  0  +1], -1:1, n, n);

			% stack matrices for coupled odes
			Z = zeros(n,n);
			M = [ M Z; Z M];
			D = [ D Z; Z D];

			obj.M  = M;
			obj.Dc = D;
		end

		% Lax-Friedrich method single time step 
		function q = advect(obj, t, q, dt)
			qold = q;
		        %f = feval(obj.pde.flux, t, q);
		        f = obj.pde.flux(t, q);
		%        p = 0.0; % choose 0 <= p < 1, if dt small
		%	 in this way, the smoothing is reduced
		%        u = p*u + (1-p)*M*u - 0.5*dt/dx*(Dc*f);
		        q = 0.5*obj.M*q - 0.5*dt/obj.dx*(obj.Dc*f);
		%	force = zeros(2*length(obj.x),1);
		%	for idx=1:length(obj.forcefun)
		%		force = force + obj.forcefun{idx}(t, obj.x, q);
		%	end
		%	q = q + dt*0.5*obj.M*force; % TODO averaging necessary here?
			%q = q + dt*force; % TODO averaging necessary here?
	        	% apply boundary condition
        		q = obj.apply_bc(t,q); %, q, qold, dt);
		end % step
	end % methods
end % class




