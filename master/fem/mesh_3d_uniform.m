% Sat Feb 25 17:54:38 MSK 2012
% Karl Kästner, Berlin

function [P T B] = mesh_3d_uniform(n, L, x0)
	if (nargin() < 3)
		x0 = [0 0 0];
	end
	% there must be at least two points per dimension
	if (any(n<2))
		warning('Increasing number of points to 2 per dimension');
		n = max(n,2);
	end

	% quick fix
	n = n([3 2 1]);
	L = L([3 2 1]);
	x0 = x0([3 2 1]);

	% minimal meshes
%	if (-1 == n(1))
%		P = [0 0 0;
%                     1 0 0;
%                     0 1 0;
%                     0 0 1];
%		T = [1 2 3 4];
%		B = [1 2 3;
%                     1 2 4;
%                     1 3 4;
%                     2 3 4];
%		return;
%	end
%	if (-2 == n(1))
%		P = [0 0 0;
%                     1 0 0;
%                     0 1 0;
%                     0 0 1;
%                    -1 0 0];
%		P(:,1) = 0.5*P(:,1)+0.5;
%		T = [2 4 1 3;
%		     1 5 3 4];
%		B = [1 2 3;
%                     1 2 4;
%                     2 3 4;
%		     1 5 3;
%                     1 5 4;
%                     3 5 4];
%		return;
%	end

	% point coordinates per axis
	X1 = L(1)*(0:n(1)-1)'/(n(1)-1);
	X2 = L(2)*(0:n(2)-1)'/(n(2)-1);
	X3 = L(3)*(0:n(3)-1)'/(n(3)-1);

	% translate
	X1 = X1 + x0(1);
	X2 = X2 + x0(2);
	X3 = X3 + x0(3);

	% vertices
	I1 = ones(n(1),1);
	I2 = ones(n(2),1);
	I3 = ones(n(3),1);
	P = [ kron(kron(X3,I2),I1), kron(kron(I3,X2),I1), kron(kron(I3,I2),X1) ];

	% start indices of hexaeders,
	% there is one less hexaeder than points per dimension
	id1 = (0:n(1)-2)';
	I1  = ones(n(1)-1,1);
	id2 = (0:n(2)-2)';
	I2  = ones(n(2)-1,1);
	id3 = (0:n(3)-2)';
	I3  = ones(n(3)-1,1);

	% tetrahedra
	% each hexaeder has has 5 children
	% T = zeros(5*prod(n-1), 4);

	% tesselation of one hexaeder
	T0 = [   1,      2, n(1)+1, n(1)*n(2)+1;
                 2, n(1)+2, n(1)+1, n(1)*n(2)+n(1)+2;
                 2, n(1)*n(2)+n(1)+2, n(1)*n(2)+1,n(1)*n(2)+2;
              n(1)+1, n(1)*n(2)+1, n(1)*n(2)+n(1)+2, n(1)*n(2)+n(1)+1;
	         2, n(1)+2, n(1)*n(2)+1, n(1)*n(2)+n(1)+2
		];
	I0   = ones(size(T0));

	if (1)
	% tesselation of the domain into cubes
	Tc =   kron(kron(id1,I2),I3) ...
	     + kron(kron(I1,id2*n(1)),I3) ...
	     + kron(kron(I1,I2),id3*n(1)*n(2));
	% tesselation of cubes into tetras
	T = kron(Tc,I0) + kron(kron(kron(I1,I2),I3),T0);
	else
		T =   kron(kron(kron(id1,I2),I3),I0) ...
		    + kron(kron(kron(I1,id2*n(1)),I3),I0) ...
		    + kron(kron(kron(I1,I2),id3*n(1)*n(2)),I0) ...
		    + kron(kron(kron(I1,I2),I3),T0);
	end

	% tesselation of (n1-1)*(n2-1)*(n3-1) hexaeder
	% there are in total (n(1)-1)*(n(2)-1)*(n(3)-1) hexaeders
	% and 5 tetrahedra per hexaeder
end % triangulate_3d

