% Jan 29 20:41 2011 (MST)
% Karl Kästner, Berlin

% todo - select not to return q

% core lanczos routinge computes the hessenberg factorisation
% Q'AQ = T, with T triagonal, Q orthogonal
% todo, parameter Q and T for continuing
function [Q T q_mp1 beta_mp1 alpha beta] = lanczos(A, m, Q0, fullq, abstol)
	if (nargin() < 5 || isempty(abstol))
		abstol = sqrt(eps);
	end
	if (nargin() < 4 || isempty(fullq))
		fullq = 1;
	end

	n = size(A,1);

	% allocate memory
	if (fullq)
		Q     = zeros(n,m+1);
	else
		Q     = zeros(n,3);
	end
	alpha = zeros(m,1);
	beta  = zeros(m+1,1);

	% initial vector
	if (isempty(Q0))
		Q(:,1) = rand(n,1)-0.5;
	else
		Q(:,1) = Q0;
	end

	% start Lancozs iteration
	beta(1) = norm(Q(:,1));
	im1 = 2; i1 = 3; i2 = 1;
	for idx=1:m
		if (fullq)
			im1 = idx-1;
			i1  = idx;
		        i2  = idx+1;
		else
			help = im1;
			im1 = i1;
			i1  = i2;
               		i2  = help;
		end

		% norm vector to 1, to avoid overflow
		Q(:,i1) = (1/beta(idx))*Q(:,i1);

		% mv-product
		if (isa(A, 'function_handle'))
			Q(:,i2) = feval(A,Q(:,i1));
		else
			Q(:,i2) = A*Q(:,i1);
		end

		% orthogonalise to first vector
		if (idx > 1)
			Q(:,i2) = Q(:,i2) - beta(idx)*Q(:,im1);
		end

		% main diagonal element - linear dependency with previous column
		alpha(idx) = Q(:,i1)'*Q(:,i2);

		% orthogonalise to second vector, remove linear depency from vector
		Q(:,i2) = Q(:,i2) - alpha(idx)*Q(:,i1);

		% off diagonal element - norm of vector
		% sqrt(q'q) is twice as fast than "norm"
		beta(idx+1) = sqrt(Q(:,i2)'*Q(:,i2));

		% stop at break down
		if (abs(beta(idx+1)) < abstol)
			dips('error, break down not yet handled');
			break
		end % if breakdown

		% recover orthogonality
		% TODO
	end % for idx
	
	% construct projected matrix
	T =  diag(sparse(beta(2:end-1)),-1) ...
	   + diag(sparse(alpha)) ...
	   + diag(sparse(beta(2:end-1)),+1);
	% slice off the last vector
	if (fullq) % && size(Q,2) > m)
		q_mp1 = Q(:,m+1)/beta(end);
		beta_mp1 = beta(end);
		Q = Q(:,1:m);
	else
		beta_mp1 = beta(end);
		q_mp1 = [];
	end
end % lanczos

