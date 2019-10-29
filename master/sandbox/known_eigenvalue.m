% Fri Sep 16 02:29:52 MSD 2011
% Karl Kästner

% for knwon eigenvector simply compute lambda = (A(1,:)*x)./x(1) is in (O(n))
% run-time: O(n^2)
function x = known_eigenvalue(A, lambda)
	% reformulate eigenvalue problem Ax = lambda x <=> (A - lambda*I)x = 0
	A = A - lambda*eye(size(A));
	% choose first coordinate to be 1
	x = zeros(length(A),1);
	x(1) = -1;
	% deflate coordinate system by 1 equation
	% and solve for remaining coordinates
	x(2:end,1) = A(2:end,2:end) \ A(2:end,1);
	% normalise
	norminv = 1/norm(x);
	x = x*norminv;
end
 
