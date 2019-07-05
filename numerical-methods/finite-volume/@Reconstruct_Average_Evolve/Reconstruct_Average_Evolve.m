% 2011-02-20 15:55:20.000000000 +0100
% \date Sun Mar 13 01:25:56 MSK 2011
% 2013/04/14 23:45
% Sa 6. Feb 23:18:07 CET 2016
% Karl Kästner, KTH Stockholm, Berlin
%% Reconstruct Average Evolve Finite Volume Method for treatment of 1+1D pdes
%%
%% McCronack Scheme
%% err = O(dt^2) + O(dx^2), except as discontinuities
%% error:
%%	h_xxx(3:end-2) = 1/dx^3*( -0.5*h(1:end-4) + h(2:end-3) - h(4:end-1)  + 0.5*h(5:end) );
%%	th = -1/6*dx^2*qh_.*(1 - (qh_*dt/dx).^2).*h_xxx;
classdef Reconstruct_Average_Evolve < Finite_Volume
%classdef Lax_Wendroff < Finite_Volume
	properties
		% averaging matrices
		Mc
		Ml
		Mr
		% difference matrices
		Dl
		Dr
		Dc
		% shift matrices
		Right
		Left
		% stepper
		advect
		limiter
		%step = @(obj,t,q,dt) Reconstruct_Average_Evolve.step_highres(obj,t,q,dt);
	end % properties
	methods
		function obj = Reconstruct_Average_Evolve()
			obj.advect = @obj.advect_highres;
		end

		function obj = init(obj,X,n)
			obj = init@Finite_Volume(obj, X,n);
			m   = obj.pde.m;
			L   = diff(X);
			obj.dx = L/(n-1);
			% TODO end boundaries of difference matrices

			% averaging matrices, TODO make convention of 0.5 consistend
			%Mc = spdiags(ones(n,1)*[0.5  0 0.5], -1:1, n, n);
			%Ml = spdiags(ones(2*n,1)*[ 1 1], -1:0, n, n);
			%Mr = spdiags(ones(2*n,1)*[ 1 1],  0:1, n, n);

			linear = true;
			% Difference Matrices
		
			if (linear)
%			Dc = spdiags(ones(n,1)*[-1, 0,+1], -1:1, n, n);
%			Dl = spdiags(ones(n,1)*[ -1 +1   0], -1:1, n, n);
			Dl = spdiags([[-1,-1, 0]; [-1,+1,+1]; ones(n-2,1)*[-1,+1, 0]], -1:1, n, n);
			%Dr = spdiags(ones(n,1)*[  0 -1  +1], -1:1, n, n);
			Dr = spdiags([ones(n-2,1)*[ 0,-1,+1]; [-1,-1,+1]; [ 0,+1,+1]], -1:1, n, n);
			%D2 = spdiags(ones(n,1)*[1 -2 1], -1:1, n, n);

			% shift matrices
			Left  = spdiags([[1,2,0];[1,0,-1]; ones(n-1,1)*[1,0,0]],-1:1,n,n);
			%Left  = spdiags([[1 1]; ones(n-1,1)*[1,0]],-1:0,n,n);
			Right = spdiags([ones(n-2,1)*[0,0,1]; [-1,0,1]; [0,2,1]],-1:1,n,n);
			%Right = spdiags([ones(n-1,1)*[0,1]; [1 1]],0:1,n,n);
			else
				Dl = spdiags([[-1,0, 0]; [-1,+1,0]; ones(n-2,1)*[-1,+1, 0]], -1:1, n, n);
				Dr = spdiags([ones(n-2,1)*[ 0,-1,+1]; [-1,-1,+1]; [ 0,+1,+1]], -1:1, n, n);

				% shift matrices
				Left = spdiags([[1,1,0];[1,0,0]; ones(n-1,1)*[1,0,0]],-1:1,n,n)
				Right = spdiags([ones(n-2,1)*[0,0,1]; [0,0,1]; [0,1,1]],-1:1,n,n);
			end

			% stack matrices for systems
			%Z = zeros(n,n);
			I         = speye(m);
%			obj.Dc    = kron(I,Dc)
			obj.Dl    = kron(I,Dl)
			obj.Dr    = kron(I,Dr)
%			obj.Mc    = kron(I,Mc);
%			obj.Ml    = kron(I,Ml);
%			obj.Mr    = kron(I,Mr);
			obj.Right = kron(I,Right);
			obj.Left = kron(I,Left);
		end % init
	end % methods
end % class RAE

