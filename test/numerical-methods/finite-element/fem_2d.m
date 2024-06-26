% Thu Feb 23 01:55:50 MSK 2012
% Karl Kästner, Berlin
%
%
% P : points
% T : triangles
% int : 2d integration scheme
% fun : potential function
%
function [A B] = fem_2d(P, T, Bc, int, fun)
	path(path,'potential');
	
	% assemble Laplacian matrix
	A = assemble_2d_3_dphi_dphi(P, T);
	% assemble potential matrix
	V = assemble_2d_3_phi_phi(P, T, int, @(q) feval(fun,q));
	% stiffness matrix
	A = -A+V;
	% assemble mass matrix
	B = assemble_2d_3_phi_phi(P, T, int);
% bug here
%	B = assemble_2d_phi_phi(P, T);

%	full(B)
%	full(B_)
%	sum(sum((B - B_).^2))
%	pause(1)

	% apply the boundary conditions	
	[A B] = boundary_2d(A, B, Bc);

%	n = sqrt(size(P,1));
%	idx=[1:n (n^2-n+1):n^2 1:n:n^2 n:n:n^2];
%	A((idx-1)*(n^2) + idx) = 1e12;
	
	%display_(P-L/2*ones(size(P)),T)
	%plot(P(:,1)-L/2,P(:,2)-L/2,'.r'); hold on
	%plot(P(idx,1),P(idx,2),'.y');
end % fem_2D

