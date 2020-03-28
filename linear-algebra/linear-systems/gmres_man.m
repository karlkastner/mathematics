% Sat Nov 20 14:16:02 MSK 2010
% Karl Kästner, KTH Stockholm

% Ax = b (LES)
% TODO do only one arnoldi and qr-factorisation step per iteration
function [x R T] = gmres_man(A,b)
	tol = 1e-12;
	nb = norm(b);
	idx = 1;
	R = [];
	i_max = 100;
	while (1)
		% H Q = Q A
		[Q H] = arnoldi(A,b,idx);
		% AQ(:,1:end-1) - Q*H = 0 (Arnoldi scheme)
		% minimise || H*y - ||b||e1 ||
		%       => y = H'H / H'||b||e1 (numerically unstable)
		%          x = Q*y;
		% y = (H'*H) \ (H'*(norm(b)*e1));
		e1 = [ 1; zeros(size(H,2),1)];
		y = H \ (nb*e1);
		% get solution vector
		x = Q(:,1:idx)*y;
		
		%% break on convergence
		%if ( ||H_k y_k - \beta e_1|| < tol)
		R(idx) = norm(H*y - nb*e1);
		if (R(idx) < tol || idx+1 > i_max)
			break;
		end
		idx = idx+1;
	end
end % function gmres_man


