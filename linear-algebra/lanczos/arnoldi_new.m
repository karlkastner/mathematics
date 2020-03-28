% Sat Sep 10 17:24:26 MSD 2011
% Karl Kästner

% test: 0 == norm(A*Q(1:16,1:16) - Q*H)
% AQ_n = Q_n+1*H_n

function [H Q W] = arnoldi(A, n, q0)
	tol = sqrt(eps);
	n = floor(n);
	if (issparse(A))
		H = spalloc(n,n-1,3*n);
	else
		H = zeros(n,n-1);
	end
	if (nargin()<3)
		q0     = rand(size(A,1),1);
	end

	Q     = zeros(length(A),n);
	%Q(:,1) = ones(length(A),1)/sqrt(length(A));	% does not work for D2
	%Q(1,1) = 1; % worse accuracy ?	% does not work for F
	Q(:,1) = q0;
	Q(:,1) = Q(:,1) / norm(Q(:,1));
	W = eye(n,n);
	%W = zeros(n,n);
	err = 0;
	alpha = zeros(n-1,1);
	beta = zeros(n,1);

	% start arnoldi iteration
	for idx=2:n
		Q(:,idx) = A*Q(:,idx-1);
		% remove linear dependencies with previous columns
		for jdx=1:idx-1
			% linear dependency with previous vector
			h = Q(:,jdx)'*Q(:,idx);
			H(jdx,idx-1) = h;
			% remove linear dependency from current vector
			Q(:,idx) = Q(:,idx) - h*Q(:,jdx);
		end % jdx
		% main diagonal element - norm of vetor
		h = norm(Q(:,idx));

		% stop at break down
		if (abs(h) < tol)
			return;
		end
		% norm vector to 1 to avoid overflow
		H(idx,idx-1) = h;
		Q(:,idx) = Q(:,idx) / h;

	end % for idx
end % function

