% Wed Feb  5 17:47:58 WIB 2014
% Karl Kastner, Berlin
%
% fem2D::interpolation_matrix.m
% 
% sets up the interpolation matrix
%
% val  = A*c
% val  : interpolated function values at the measured points
% c    : coefficient of the polynomial
% A    : vandermonde interpolation matrix to compute the function values at the measured points from the polynomial coefficients
% val_p: function values at the vertices
% A_p  : vandermonde interpolation matrix for computing the function values at the vertices from the polynomial coefficinets
% val_p = A_p *c
% regression to obtain values at vertices :
% val = A*c = A*(A_p^-1)*val_p = AA*val_p
% A_p (A^-1) val = val_p
%
% TODO this sets up the vertex and measurement matrix in one step, better deentangle ?
% TODO kick out ill defined elements/points for elements with insufficiently many measurement points inbetween
%
function A = interpolation_matrix(mesh,pdx,x,y)
	nbuf = 0;
	buf = zeros(3*length(x),3);
	% for each element
	for edx = 1:size(mesh.T,1)
		% grep the point coordinates
		x_p = mesh.P(mesh.T(edx,:),1);
		y_p = mesh.P(mesh.T(edx,:),2);

		% get indices of measurement points that are within this element
		fdx = pdx{edx};

		% TODO, for quads, that should be at least 4
		if (length(fdx) >= 3)

		switch (mesh.type)
		case {'triangle'}
			% set up the vertex vandermonde matrix
			Ap = [ 1 x_p(1) y_p(1)
			       1 x_p(2) y_p(2)
		               1 x_p(3) y_p(3) ];
			% set up the measurement vandermonde matrix
			A  = [ ones(length(fdx),1) x(fdx) y(fdx)];
			ndim = 3;
		case {'quadrilateral'}
			% set up the vertex vandermonde matrix
			Ap = [ 1 x_p(1) y_p(1) x_p(1)*y_p(1)
			       1 x_p(2) y_p(2) x_p(2)*y_p(2)
		               1 x_p(3) y_p(3) x_p(3)*y_p(3)
			       1 x_p(4) y_p(4) x_p(4)*y_p(4) ];
			% set up the measurement vandermonde matrix
			A  = [ ones(length(fdx),1) x(fdx) y(fdx) x(fdx).*y(fdx) ];
			ndim = 4;
		otherwise
			error('element type not yet implemented');
		end

		% combined vandermonde matrix
		% AA interpolates from values at the vertices to values at the measured points
		% TODO : precalculated inverse of the vandermonde matrix (in n^2 instead of n^3 steps)
		AA = A*inv(Ap);
		%AA = Ap \ A;

		% write into global matrix
		for idx=1:length(fdx)
		 for jdx=1:ndim
		  nbuf = nbuf+1;
		  buf = resizebuf(buf,nbuf);
		  buf(nbuf,1) = mesh.T(edx,jdx);
		  buf(nbuf,2) = fdx(idx);
		  buf(nbuf,3) = AA(idx,jdx);
		 end % for idx
		end % for jdx
		end % if length(fdx) >= 3
	end % for edx
	
	% set up regression matrix from buffer
	n = size(mesh.P,1);
	m = length(x);
	A = sparse(buf(1:nbuf,1), buf(1:nbuf,2), buf(1:nbuf,3),n,m);
end % function interpolation_matrix()

function buf = resizebuf(buf,nbuf)
	if (nbuf > size(buf,1))
		buf(2*end,1) = NaN;
	end
end

