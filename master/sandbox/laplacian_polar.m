% Laplacian in polar coordinates:
%R = diag(1./r(2:N2+1));
%Z = zeros(M2);
%I = eye(M2);
%L = kron(D1+R*E1, eye(M)) + kron(D2+R*E2, [Z I;I Z]) + kron(R^2, D2t);
% Numerical Analysis of Partial Differential Equations AvShiu-Hong Lui,S. H, Lui

% returns transformed laplacian r psi = phi
function [L theta r L_] = laplacian_polar(n,L0)
	 hr = 1/(n+1);
	 r = hr*L0*(1:n);
	 ht = 2*pi/n;
	 theta = ht*(0:n-1);
	 Ri = diag(sparse(1./r));
	 I = speye(n);

	 rr = r+0.5*r(1); %0.5*(r + [r(2:end)      (2*r(end)-r(end-1))]);
	 rl = r-0.5*r(1); %0.5*(r + [(2*r(1)-r(2)) r(1:end-1)]);
	 D2 = 1/hr^2*(- diag(sparse(rl))*spdiags(ones(n,1)*[-1 1],-1:0,n,n) ...
		      + diag(sparse(rr))*spdiags(ones(n,1)*[-1 1], 0:1,n,n) );
	
	 dR  = 0.5/hr*spdiags(ones(n,1)*[-1 0 1],-1:1,n,n);
	 ddR =   1/hr^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	 ddT =   1/ht^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	 % periodic boundary conditions for theta
	 ddT(1,end) = ddT(1,2);
	 ddT(end,1) = ddT(end,end-1);

%F1=full(Ri*D2);
%F1(1:8,1:8)
%F2=full(Ri*dR + ddR);
%F2(1:8,1:8)
%pause

	 L   = kron(D2,I)          + kron(Ri, ddT);
 	 L_  = kron(Ri*dR + ddR,I) + kron(Ri.^2,ddT);
end % laplacian_polar

