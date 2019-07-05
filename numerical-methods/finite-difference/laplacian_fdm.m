% Di 26. Jan 11:57:00 CET 2016
% Karl Kastner, Berlin
%
%% finite difference matrix of the laplacian
function A = laplacian_fdm(n,L,bc)
	A = laplacian_fdm1(n(1),L(1),bc);
	for idx=2:length(n)
		Ai = laplacian_fdm1(n(idx),L(idx),bc);
		if (idx > 1)
			I  = speye(size(A,1));
			Ii = speye(n(idx));
			A  = kron(A,Ii) + kron(I,Ai);
		end
	end
end

function A = laplacian_fdm1(n,L,bc)
	switch (n)
	case {0}
		A = [];
	case {1}
		A = 0; % 1 for dirichlet?
	otherwise
	h = L/n;
	% difference matrix
	A = (1/h^2)*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);

	% apply BC
	switch (lower(bc))
	case {'dirichlet'}
		% homogeneous Dirichlet
		% remove first and last column, because f(1) = f(n) = 0
		% (setting first row to all zeros yields singular system)
%		A(1,:)   = []; A(1,1:2) = [1,0]
%		A(end,:) = [];
	case {'neumann'}
		% homogeneous Neumann
		% left side with ghost point:
		%   1 | -2 1  = 1 % pde
		%  -1 |  1 0  = 0 % bc
		% ----------
		%   0 | -1 1  = 1 % left boundary w/o ghost point 
		A(1,1:2)         = 1/h^2*[-1 1];
		A(end,end-1:end) = 1/h^2*[1 -1];
	otherwise
		error('unimplmented Boundary condition');
	end
	end
end
%		
%	
%
%% BC
%bc = 3;
%switch (bc)
%case {0}
% A(1,:) = 0;
% A(1,2) = 1;
% A(end,:) = 0;
% A(end,end-1) = -1;
%case {1}
%	A(:,1) = []; A(1,:) = [];
%	A(:,end) = []; A(end,:) = [];
%case {2}
%	A(:,1) = []; A(1,:) = [];
% A(1,1) = -1/h^2;
%	A(:,end) = []; A(end,:) = [];
% A(end,end) = -1/h^2;
%case {3}
%	%A(1,1:3)         = -2/h^2*[1,-2,1];
%	%A(end,end-2:end) = -2/h^2*[1,-2,1];
%	%A(1,1:3)         = [1,-2,1];
%	%A(end,end-2:end) = [1,-2,1];
%	%A(1,1:3)         = [0,1,0];
%	%A(end,end-2:end) = [0,1,0];
%	A(1,1:2)         = 1/h^2*[-1 1];
%	A(end,end-1:end) = 1/h^2*[1 -1];
%	
%	% symmetry trick
%	%A(2,1:2) 	   = A(2,1:2) + ([1 -1] - [1 -1]/h^2);
%	%A(end-1,end-1:end) = A(end-1,end-1:end) + ([-1 1] - [-1 1]/h^2);
%end
%
