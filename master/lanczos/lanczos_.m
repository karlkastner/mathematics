% Sat Sep 10 17:24:26 MSD 2011
% Karl Kästner

% test: 0 == norm(A*Q(1:16,1:16) - Q*H)
% AQ_n = Q_n+1*H_n

function [Q H W] = lanczos(A, n, recover)
	if (nargin < 3)
		recover = 0;
	end
	if (nargin < 2)
		n = size(A,1);
	end

	tol = 1e-7;
	n = min(size(A,1),floor(n));

%	if (issparse(A))
%		H = spalloc(n,n-1,3*n);
%	else
%		H = zeros(n,n-1);
%	end

%	Q(:,1) = ones(length(A),1);	% does not work for D2 ?
%	Q = zeros(length(A),n); Q(1,1) = 1; % worse accuracy ?	% does not work for F
	Q(:,1) = rand(length(A),1);
	Q(:,1) = Q(:,1) / norm(Q(:,1));
	W = eye(n,n);
	%W = zeros(n,n);
	err = 0;
%	alpha = zeros(n-1,1);
%	beta = zeros(n,1);
%	q_old = zeros(size(A,1),1);

	beta = 0;
	% start Lancozs iteration
	for idx=1:n
		if (idx > 1)
			Q(:,idx+1) = A*Q(:,idx) - beta(idx)*Q(:,idx-1);
			%q_old = Q(:,idx-1);
			%Q(:,idx+1) = A*Q(:,idx) - beta(idx)*q_old;
		else
			Q(:,idx+1) = A*Q(:,idx);
		end

		% main diagonal element - linear dependency with previous column
		alpha(idx) = Q(:,idx)'*Q(:,idx+1);

		% remove linear depency from vector
		Q(:,idx+1) = Q(:,idx+1) - alpha(idx)*Q(:,idx);

		% off diagonal element - norm of vector
		beta(idx+1) = norm(Q(:,idx+1));

		% stop at break down
		if (abs(beta(idx+1)) < tol)
			break;
		end % if breakdown

		% norm vector to 1, to avoid overflow
		betainv = 1/beta(idx+1);
		Q(:,idx+1) = Q(:,idx+1) * betainv;

		% recover orthogonality
			switch (recover)
				case {0}
				% no reothogonalisation
				case {1}
				% full reorthogonalisation
					for jdx = 1:idx
						Q(:,idx+1) = Q(:,idx+1) - (Q(:,idx+1)'*Q(:,jdx))*Q(:,jdx);
					end
				case {2}
				% periodic reorthogonalisation
					[W err] = update_W(W, idx, alpha, beta, Q);
					while (err > sqrt(eps))
						[Q W beta] = reorthogonalise(Q, idx, W, beta);

						% todo - store only last two rows of W
						% todo compute err in function reo
						% todo correct beta and beta -1 and renorm vectors ??
						% 	q_p1 = b_p1 * q_p1; q_p1 = q_p1 / norm (q_p1)


						%Q(:,idx+1) = Q(:,idx+1)/norm(Q(:,idx+1));

						err = max(abs(W(idx+1,1:idx)));
					end % while
				otherwise
					msgID = 'arnoldi:lanczos';
					msg = 'unknown reorthogonalisation option';
					error(msgID, msg);
			end % switch method
	end % for idx
	% construct matrix
%	Q
	H = sparse(diag(alpha) + diag(beta(2:end-1),-1) + diag(beta(2:end-1),+1));
	%H(n,n-1) = beta(end);
	if (nargout < 2)
		Q = H;	
	end
end % lanczos

% construct W by recurrence
function [W err] = update_W(W, idx, alpha_, beta_, Q)
	% main diagonal	
        %W(idx+1,idx+1) = 1;

	% round off error of last vector
	w = Q(:,idx+1)'*Q(:,idx);
	err = w;
	W(idx+1,idx) = w;
	W(idx,idx+1) = w;

	% first column (propagated round off error)
	if ( idx+1 > 2 )
	     w = 1/beta_(idx+1) * (   beta_(2) * W(idx,1+1) ...
				   + (alpha_(1) - alpha_(idx))*W(idx,1) ...
                                   - beta_(idx)*W(idx-1,1) );
             W(idx+1,1) = w;
	     W(1,idx+1) = w;
	     if (abs(w) > err)
		err = abs(w);
	     end
	end % first column

	% remaining columns
	for kdx=2:idx-1
		w = 1/beta_(idx+1) * ( beta_(kdx+1) * W(idx,kdx+1) ...
					   + (alpha_(kdx) - alpha_(idx))*W(idx,kdx) ...
                                           + beta_(kdx)*W(idx,kdx-1) ...
                                           - beta_(idx)*W(idx-1,kdx) );
                W(idx+1,kdx) = w;
		W(kdx,idx+1) = w;
		if (abs(w) > err)
			err = abs(w);
		end
        end % remaining columns
	%err = max(abs(W(1:idx,idx+1)));
end % update_W

% restores orthogonlatity of columns idx and idx+1
function [Q W beta] = reorthogonalise(Q, idx, W, beta)
	% reortogonalise last two vectors
	for jdx = 1:idx-1
		% using values in w instead of recalculation Q'Q is better !!!
		%Q(:,idx) = Q(:,idx) - (Q(:,idx)'*Q(:,jdx))*Q(:,jdx);
		%Q(:,idx+1) = Q(:,idx+1) - (Q(:,idx+1)'*Q(:,jdx))*Q(:,jdx);
		Q(:,idx) = Q(:,idx) - W(idx, jdx)*Q(:,jdx);
		Q(:,idx+1) = Q(:,idx+1) - W(idx+1, jdx)*Q(:,jdx);
	end
	%Q(:,idx+1) = Q(:,idx+1) - (Q(:,idx+1)'*Q(:,idx))*Q(:,idx);
	Q(:,idx+1) = Q(:,idx+1) - W(idx+1,idx)*Q(:,idx);
	
	% renorm last two vectors - worse results !!!
%{
	 norminv = 1/(Q(:,idx)'*Q(:,idx));
	 Q(:,idx) = Q(:,idx)*norminv;
	 norminv = 1/(Q(:,idx+1)'*Q(:,idx+1));
	 Q(:,idx+1) = Q(:,idx+1)*norminv;
%}
%{ 
	Q(:,idx)  = beta(idx)*Q(:,idx);
	beta(idx) = Q(:,idx)'*Q(:,idx);
	Q(:,idx)  = Q(:,idx)/beta(idx);

	if(idx+1 > 2)
		Q(:,idx+1)  = beta(idx+1)*Q(:,idx+1);
		beta(idx+1) = Q(:,idx+1)'*Q(:,idx+1);
		Q(:,idx+1)  = Q(:,idx+1)/beta(idx+1);
	end
%}
	
	err = 0.0;

	% recalculate dependencies
	for jdx = 1:idx-1
		w = Q(:,idx)'*Q(:,jdx);
		W(idx,jdx) = w;
		W(jdx,idx) = w;
		if (abs(err) > w)
			err = w;
		end
		w = Q(:,idx+1)'*Q(:,jdx);
		W(idx+1,jdx) = w;
		W(jdx,idx+1) = w;
		if (abs(err) > w)
			err = w;
		end
	end
	w = Q(:,idx+1)'*Q(:,idx);
	W(idx+1,idx) = w;
	W(idx, idx+1) = w;
	if (abs(err) > w)
		err = w;
	end

end % reorthogonalis

