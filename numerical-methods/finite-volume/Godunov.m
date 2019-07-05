% Thu Apr 28 13:00:19 MSD 2011
% Thu Apr 28 14:15:20 MSD 2011
% Sa 6. Feb 23:25:59 CET 2016

%% Godunov, upwind method for systems of pdes
classdef Godunov < Finite_Volume
	properties
		fluxmat
%		x
%		dx	

		% difference matrices
		Dl
		Dr
		bc
		hrflag
	end % properties
	methods
	function obj = init(obj, X, n)
		obj = init@Finite_Volume(obj, X,n);
	
		Z  = spalloc(n,n,0);
		Dl = spdiags(ones(n,1)*[-1 +1  0], -1:1, n, n);
		Dr = spdiags(ones(n,1)*[ 0 -1 +1], -1:1, n, n);

		obj.Dr = [Dr Z; Z Dr];
		obj.Dl = [Dl Z; Z Dl];
	end
	function [q obj] = step(obj, t, q)

		% get eigenvalues
		[R Lambda Rinv] = feval(obj.fluxmat, t, q);

		% 12.7
		w_l = obj.Dl*q;
		w_r = obj.Dr*q;
		s_l = obj.Dl*f./w_l;	% what if 0 == w_l ??
		s_r = obj.Dr*f./w_l;	% what if 0 == w_l ??
		s = (f(qr) - f(ql)) / (qr - ql); % 11.21

		% shure ?
		sp_l = (s_l > 0) .* s_l;
		sn_r = (s_r < 0) .* s_r;

		% 12.8 - alternative see 15.6
		% see 15.8, 15.10
		% Ahat_l = f'(qhat_l) = f'(obj.Dl*q) % near 15.15
		% Ahat_l = 0.5*(f'(ql) + f'(q)]
		ApdQl = sp_l.*w_l;
		AndQr = sn_r.*w_r;

		% 12.5, 15.5
		q = q - dt/dx(ApdQl - AndQr);
	
		% split into left and right eigenvalues
		Lambda_m = (Lambda < 0).*Lambda;
		Lambda_p = (Lambda > 0).*Lambda;

		% left and right matrix
		A_l = R*Lambda_m*Rinv;
		A_r = R*Lambda_p*Rinv;
	
		if (1)
			% F = Ar ql + Al q
			% q = q + dt/dx Dl F

			% variant 1 : 4.43, 4.47
			q = q - obj.dt/obj.dx*(A_r*(obj.Dl*q) + A_l*(obj.Dr*q));
		else
			% variant 2 : 4.51 (ident), 4.56
			w_l   = R*diag(Rinv*obj.Dl*sparse(q));
			w_r   = R*diag(Rinv*obj.Dl*sparse(q));
			q = q - dt/obj.dx*(w_r.*Lambda_m + w_l.*Lambda_p);
		end
	
		% apply boundary condition	
		q = feval(bc, t, q);
	end % step

	end % methods
end % godunov

