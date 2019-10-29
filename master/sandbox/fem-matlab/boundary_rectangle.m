function [q,g,h,r]=boundary(p,e,u,time)

	length(e)
	% Neumann
	%q = zeros(1,2*length(e));
	%g = zeros(1,2*length(e));
	N = 1; %size(p,2);
	ne = size(e,2)
	
	q = zeros( N^2, size(e,2));
	g = zeros(   N, size(e,2));

	% Dirichlet hu=r
	h =  ones(N^2,2*ne);
	r = zeros(N,2*ne);

end % boundary_rectangle

