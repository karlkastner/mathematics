% 2012 Feb 23 02:17
% Karl Kästner, Berlin
%
% triangulate a two dimensional square domain
% n+1 : number of elements pre row and column
% L0 : width and height of domain

function [P T B X] = triangulation_uniform(arg1, arg2, x0)
	if (nargin() < 3)
		x0 = [0 0];
	end
	if (isempty(arg1))
		% arg1 is ignored
		% arg2 is cell array with X and Y vectors
		X1 = arg2{1};
		X2 = arg2{2};
		X = {X1, X2};
		n(1) = length(X1);
		n(2) = length(X2);
	else
		% arg1(1,2) is number of points per dimension
		n = arg1;
		n = max(n,2);
%		if (-1 == n(1))
%			P = [0 0;
%	                     1 0;
%	                     0 1];
%			T = [1 2 3];
%			B = [  1 2 100;
%	                       2 3 200;
%	                       3 1 300]; 
%			return
%		end
		% arg2 is length of domain per dimension
		L0 = arg2;
		h  = L0./(n-1);
	end
	L0 = [L0(2),L0(1)];
	n  = [n(2),n(1)];
	x0 = [x0(2),x0(1)];

	% point coordinates per axis
	X1 = L0(1)*(0:n(1)-1)'/(n(1)-1);
	X2 = L0(2)*(0:n(2)-1)'/(n(2)-1);

	% translate
	X1 = X1 + x0(1);
	X2 = X2 + x0(2);

	% vertices
	I1 = ones(n(1),1);
	I2 = ones(n(2),1);
	P = [ kron(X2,I1), kron(I2,X1) ];

	% start indices of rectangles,
	% there is one less rectangle than vertices per dimension
	id1 = (0:n(1)-2)';
	I1  = ones(n(1)-1,1);
	id2 = (0:n(2)-2)';
	I2  = ones(n(2)-1,1);

	% triangles
	% each rectangle has 2 children

	% tesselation of one rectangle
	T0 = [   1,      2, n(1)+1;
                 2, n(1)+2, n(1)+1];
	I0   = ones(size(T0));

	% tesselation of the domain into rectangles
	Tc =   kron(id1,I2) ...
	     + kron(I1,id2*n(1));

	% tesselation of rectangles into triangles
	T = kron(Tc,I0) + kron(kron(I1,I2),T0);

	B = [];
end % mesh_2d_uniform



