% Mon 13 Jun 09:39:34 CEST 2022
% A^p = V l^p V'
%
% TODO allow for different boundary conditions
%
% to compute the laplacian matrix : x = eye(prod(n))
function z = laplacian_power(L,n,p,x)
	y = zeros(size(x));

	% right vectors and powers
	% note that right vectors can be skipped for random x,y
%	x = reshape(x,n);
	k = 1;
	for i=1:n(1)
	 for j=1:n(2)
		[vx,vy] = laplacian_eigenvector(n,[i,j]);
%		v = kron(vy,vx);
		v = vy*vx';
		v = flat(v);
		y(k,:) = v'*x;
		k = k+1;
	 end
	end

	% diagonal
	k = 1;
	for i=1:n(1)
	 for j=1:n(2)
		e = sum(laplacian_eigenvalue(L,n,[i,j]));
		y(k,:) = (e).^p*y(k,:);
		k = k+1;
	 end
	end

	% left vectors
	z = zeros(size(x));
	k = 1;
	for i=1:n(1)
	 for j=1:n(2)
		[vx,vy] = laplacian_eigenvector(n,[i,j]);
		v = vy*vx';
%		v = v';
		v = flat(v);
		%v = flat(v);
		%v = kron(vx,vy);
		z(k,:) = v'*y;
		k = k+1;
	 end
	end
end

