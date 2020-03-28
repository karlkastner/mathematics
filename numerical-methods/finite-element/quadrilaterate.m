% Tue Jan 28 15:09:30 WIB 2014
% Karl Kastner, Berlin

% meshes the unit square into equal quadrilaterals
%
% n(1): number of points on x axis
% n(2): number of points on y axis
%
% Q : rows define quadrilateral corners as indices in point matrix P
% P : point coordinates
function [Q P] = quadrilaterate(n)
	X = linspace(0,1,n(1))';
	Y = linspace(0,1,n(2))';
	P = [kron(ones(n(2),1),X) kron(Y,ones(n(1),1))];
	Q = [];
	Q_ = [(1:n(1)-1)' (2:n(1))' (n(1)+2:2*n(1))' (n(1)+1:2*n(1)-1)'];
	for jdx=1:n(2)-1
		Q = [Q; Q_ + (jdx-1)*n(1)];
	end
%	Q = Q'; <- for patch, but not to be consistend with the mesh class!!!
end

