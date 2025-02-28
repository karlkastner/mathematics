% Wed  6 Nov 11:16:03 CET 2024
function [dy_de,D2] = grad_diffuse_1d_fdm(y,dt,e,n,L)
	D2 = derivative_matrix_2_1d(n,L,2,'circular');
	I  = speye(n);
	A  = (I-e*D2);
	dy_de = ((A*A)\(D2*y));
	% equivalent to:
	% dy_de = D2*((A*A) \ y);
end

