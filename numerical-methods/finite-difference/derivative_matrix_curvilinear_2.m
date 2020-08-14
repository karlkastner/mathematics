% Wed 29 Aug 11:03:03 CEST 2018
%
%% derivative matrix on a two dimensional curvilinear grid
%% the grid has not necessarily to be orthogonal
% TODO mixed derivative is for some reason buggy
function [Dx,Dy,D2x,Dxy,D2y,L] = derivative_matrix_curvilinear_2(x,y,isorthogonal)
	if (nargin()<3)
		isorthogonal = false;
	end

	n = size(x);
	x = flat(x);
	y = flat(y);

	% unit difference matrices
	[Ds,Dn,D2s,Dsn,D2n] = derivative_matrix_2d(n,n-1);

	% forward transformation matrix and step width
	x_s = Ds*x;
	x_n = Dn*x;
	y_s = Ds*y;
	y_n = Dn*y;

	x_ss = D2s*x;
	x_sn = Dsn*x;
	x_nn = D2n*x;
	y_ss = D2s*y;
	y_sn = Dsn*y;
	y_nn = D2n*y;

	J = x_s.*y_n - x_n.*y_s;
	a = (x_n.*y_sn + x_nn.*y_s - x_s.*y_nn - x_sn.*y_n)./J.^2;
	b = (x_n.*y_ss - x_ss.*y_n - x_s.*y_sn + x_sn.*y_s)./J.^2;

	Dx = diag(sparse(-y_s./J))*Dn + diag(sparse( y_n./J))*Ds;
	Dy = diag(sparse( x_s./J))*Dn + diag(sparse(-x_n./J))*Ds;

	D2x = ( ...
	  diag(sparse( (y_n.^2./J.^2)))*D2s ...
	- diag(sparse((2.*y_n.*y_s./J.^2)))*Dsn ...
	+ diag(sparse((y_s.^2./J.^2)))*D2n ...
	- diag(sparse((y_s.*(y_nn./J + y_n.*a) - y_n.*(y_sn./J + y_n.*b))./J))*Ds ...
	+ diag(sparse((y_s.*(y_sn./J + y_s.*a) - y_n.*(y_ss./J + y_s.*b))./J))*Dn ...
	);
	Dxy = ( ...
	-diag(sparse((x_n.*y_n)./J.^2))*D2s ...
	+diag(sparse((((x_n.*y_s)./J + (x_s.*y_n)./J))./J))*Dsn ...
	-diag(sparse((x_s.*y_s)./J.^2))*D2n ...
	+diag(sparse(((y_s.*(x_nn./J + (x_n.*a)) - y_n.*(x_sn./J + (x_n.*b))))./J))*Ds ...
	-diag(sparse(((y_s.*(x_sn./J + (x_s.*a)) - y_n.*(x_ss./J + (x_s.*b))))./J))*Dn ...
	);
	
	D2y = ( ...
	 diag(sparse( (x_s.^2)./J.^2))*D2n ...
	-diag(sparse((2.*x_n.*x_s)./J.^2))*Dsn ...
	+diag(sparse((x_n.^2)./J.^2))*D2s ...
	-diag(sparse(((x_s.*(x_nn./J + x_n.*a) - x_n.*(x_sn./J + x_n.*b)))./J))*Ds  ...
	+diag(sparse(((x_s.*(x_sn./J + x_s.*a) - x_n.*(x_ss./J + x_s.*b)))./J))*Dn  ...
	);
	
	L = D2x+D2y;

end

