% Sun 20 May 16:48:42 CEST 2018
%
% smoothing matrix for 2D grid
%   0  1/8  0
%  1/8 1/2 1/8
%   0  1/8  0 
function A = averaging_matrix_2(n)
	p = 0.5;

	A1 = spdiags(ones(n(1),1)*[1,0,1], -1:1, n(1), n(1));
	A1(1,:) = 0;
	A1(1,1:2) = [1,1];
	A1(end,:) = 0;
	A1(end,end-1:end) = [1,1];

	A2        = spdiags(ones(n(2),1)*[1,0,1], -1:1, n(2), n(2));
	A2(1,:)   = 0;
	A2(1,1:2) = [1,1];
	A2(end,:) = 0;
	A2(end,end-1:end) = [1,1];
	full(A2)

	A = (1-p)*speye(n(1)*n(2)) + p/4*(kron(A1,speye(n(2))) + kron(speye(n(1)),A2));
end

