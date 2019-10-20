% Thu Apr 28 04:30:48 MSD 2011
% Sat Apr 30 04:04:09 MSD 2011
% So 7. Feb 22:21:46 CET 2016
% So 7. Feb 13:40:12 CET 2016
% Karl Kästner

%% non linear roe solver for the SWE (randall, leveque 15.3.1)
%%
%% The roe solver guarantess:
%% - A is diagonalisable with real eigenvalues (15.12)
%% - can be determined by a closed formula
%% - is an efficient replacement for true Rieman solver
%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
classdef RoeLow < Finite_Volume
	properties
		fluxmateig
		roeaverage
		x
		dx	

		% averaging matrix
		Ml
		% difference matrices
		Dr
		Dl

		bc
		% hrflag
		limiter = @Finite_Volume.minmod;
	end % properties
	methods
		function obj = init(obj, L, n)
			obj.dx = L/(n-1);

			Z  = spalloc(n,n,0);
			Ml = 0.5*spdiags(ones(n,1)*[1 1 0], -1:1, n, n);
			Mr = 0.5*spdiags(ones(n,1)*[0 1 1], -1:1, n, n);
			Dl = spdiags(ones(n,1)*[-1 +1  0], -1:1, n, n);
			Dr = spdiags(ones(n,1)*[ 0 -1 +1], -1:1, n, n);
	
			% stack matrices for coupled odes
			obj.Dl = [Dl Z; Z Dl];
			obj.Dr = [Dr Z; Z Dr];
			obj.Ml = [Ml Z; Z Ml];
			% Mr = [ Mr Z; Z Mr];
		end
		% Roe-solver
		%TODO entropy fix
		%L = diag(L); L = diag(L - (L<0)*1E7 + (L>=0)*1E7);
		% Fri Apr 29 00:36:21 MSD 2011
		function [q obj] = step(obj, t, q)

			% low-resolution part

		        % eigenvalues of the underlying ode
			% flux = (R L Rinv)q = A q
%		[Lambda R Rinv] = feval(obj.fluxmateig, t, q, dt);

			% step upstream for positive and downstream for negative eigenvalues
%		% below 15.64
%		s_l = obj.Ml*Lambda;
%		s_r = obj.Mr*Lambda;

		        % roe averages at volume interfaces
			ql_hat = obj.roeaverage(q);
			% PDE flux matrix at volume interfaces (15.34)
			[R_hat L_hat Rinv_hat] = feval(obj.fluxmat, t, ql_hat);
			%absAl = R_l*abs(L_l)*Rinv_l;
			A_l = R_hat*abs(L_hat.*(L_hat<0))*Rinv_hat;
			A_r = R_hat*abs(L_hat.*(L_hat>0))*Rinv_hat;

			absSl = abs(L_hat);
			absSr = abs([L_hat(2:end);0]);

			% (15.64) 6.55 6.50, 4.47
			% decouplig eigenspaces
			al = Rinv_hat*obj.Dl*q;
			ar = Rinv_hat*obj.Dr*q;

			% 6.51
			al_     = left(al).*(L_hat >= 0) + right(al).*(L_hat < 0); 
			ar_     = left(ar).*(L_hat >= 0) + right(al).*(L_hat < 0); 

			% 6.51
			thetal  = al_./al;
			thetar  = ar_./ar;

			% 6.50
			% apply flux limter
			atildel = al.*obj.limiter(thetal);
			atilder = ar.*obj.limiter(thetar);

			% 6.55
			% retransformation to euclidean basis
			%w_tilde_l = minmod(w_l, w_ll*(s_l > 0) + w_r*(s_l < 0) );
			%w_tilde_r = minmod(w_r, w_l*(s_r > 0) + w_rr*(s_r < 0) );
			wtilde_l = R_hat*atilde_l;
			wtilde_l = R_hat*atilde_r;

			% flux
			% 15.63, 15.40, 6.60, 5.8, 4.61
			% note : Fr(n) = Fl(n+1) => Fr - Fl = Fl(n+1) - Fl(n) = Dr*Fl
			Fl = 0.5*(absSl.*(1 - dt/dx*absSl)*wtilde_l);
			Fr = 0.5*(absSr.*(1 - dt/dx*absSr)*wtilde_r);
			% 6.53 flux correction term
			Ftilde = 0.5*matvec2x2(,aR-dt/dx*matvec2x2(Aa,aR));
			% 6.52 flux
			F = ApQl + AmQc + Ftilde;
			
			% Fl = obj.Ml*f - 0.5*absAl*(obj.Dl*q);
			% F_tilde_l = w_tilde_l*diag((abs(s_l).*(obj.I - dt/dx*abs(s_l))));
			% F_tilde_r = w_tilde_r*diag((abs(s_r).*(obj.I - dt/dx*abs(s_r))));
	
			% here f are the fluxes at points, not at interfaces?
			% F = 0.5*(f(:,2:end) + f(:,1:end-1)) - 0.5*(Aa Dl q);

			% forcing is applied by splitting
			%source = obj.source(t,x,q);

			% finite volume update
			% 15.62, 6.59, 6.59, 5.7, 4.4
			% really Dr and Dl ???
			% q = q - dt/obj.dx*(obj.Dr*Fl);
			% high-resolution correction
			%q = q_low + 0.5*dt/dx*(F_tilde_r - F_tilde_l);
			q = q - dt/dx*(Al*(obj.Dr*q) - Ar*(obj.Dl*q)) ...
			      - dt/dx*(Fr - Fl);
                        %     + dt*source;

			% apply boundary conditions
		        q = feval(obj.bc, t, q);
	
			function x=left(x)
				%x=[x(1);x(1:end-1)];
				x = circshift(x,-1);
			end
			function x=right(x)
				%x=[x(2:end);x(end)];
				x = circshift(x,+1);
			end
		end % step
end % roeLow

