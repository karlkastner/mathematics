
function [dy_de,D2] = grad_diffuse_2d_fdm(y,dt,e,n,L)
	[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,L,2,{'circular','circular'});
	D2 = D2x+D2y;
	I  = speye(prod(n));
	A  = (I-e*D2);
	dy_de = ((A*A)\(D2*y));
	% equivalent to:
	% dy_de = D2*((A*A) \ y);
end


