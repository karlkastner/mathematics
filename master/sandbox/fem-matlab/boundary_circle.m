function [q,g,h,r]=boundary(p,e,u,time)
	%Boundary condition data
	% Neumann-conditions, zero values for Dirichlet-Boundary
	q=zeros(1,length(e));
	g=zeros(1,length(e));
	% Dirichlet-Boundary: hu=r
	h=ones(1,2*length(e));
	r=zeros(1,2*length(e));
end

