% 2012 Feb 23 02:17
% Karl Kästner, Berlin
%
% triangulate a two dimensional square domain
% n+1 : number of elements pre row and column
% L0 : width and height of domain

function [P T B X] = mesh_2d_uniform(arg1, arg2, x0)
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
		% arg1 is number of points per dimension
		n  = arg1;
	if (-1 == n(1))
		P = [0 0;
                     1 0;
                     0 1];
		T = [1 2 3];
		B = [  1 2 100;
                       2 3 200;
                       3 1 300]; 
		return
	end
		% arg2 is length of domain per dimension
		L0 = arg2;
		X1 = L0(1)/(n(1)-1)*(0:n(1)-1)';
		X2 = L0(2)/(n(2)-1)*(0:n(2)-1)';
		X = {X1, X2};
	end

	% shift
	X1 = X1 - x0(1);
	X2 = X2 - x0(2);

	I1 = ones(n(1),1);
	I2 = ones(n(2),1);
	P = [kron(X1,I2) kron(I1,X2)];
	% shift
	%P = P - ones(size(P,1),1)*x0;

	T = [];
	for idx=1:n(1)-1
		N = (idx-1)*n(2) + (1:n(2)-1)';
		T = [T;
		     N N+1   N+n(2)+1;
		     N N+n(2)+1 N+n(2)];
	end

	% x y boundary
	B = [	(1:n(2)-1)'                (2:n(2))'                      ; %    ones(n(2)-1,1);  % left
		(1:n(2):n(2)*(n(1)-2)+1)'  (n(2)+1:n(2):n(2)*(n(1)-1)+1)' ; % 2*ones(n(1)-1,1);  % bottom
		(n(2):n(2):n(2)*(n(1)-1))' (2*n(2):n(2):n(1)*n(2))'       ; % 3*ones(n(1)-1,1);  % top
		n(2)*(n(1)-1)+(1:n(2)-1)'  n(2)*(n(1)-1)+(2:n(2))'        ]; % 4*ones(n(2)-1,1)]; % right
end % mesh_2d_uniform


