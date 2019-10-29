% Sat Sep 10 17:24:26 MSD 2011
% Karl Kästner

% todo - get start vector, do not do eigs on the full system
% todo - do not calculate eigenvectors at each iteration
% todo - remove superficial alpha, beta vectors
% todo - e_true calculation out of the loop
% todo - continue the factorisation if insufficient
function [E E1 E2] = eigs_lanczos(A, k, m, sigma)
	n = size(A,1);

	if (nargin < 2)
		k = 6;
	end

	if (nargin < 3)
		sigma = 'SM';
	end

	% allocate memory
	R  = zeros(m,1);
	K  = zeros(m,1);

	% for test purposes
	if (nargout > 5)
		E_true = sort(eig(A));
		err_E = zeros(m,1);
		err_Q = zeros(m,1);
	end
	
	Q0 = ones(n,1) + 0.1*rand(n,1);
	
	% matrix factorisation
	[Q T_ void alpha beta] = lanczos(A,m,Q0,0,[]);

	% get the eigenvalues of T1
	E1 = eig_symmetric(T_,1);

	% get the eigenvalues of T2
	E2 = eig_symmetric(T_(2:end,2:end),1);

	% todo - sort out the spurious eigenvalues
	E = E1;
if (0)
	switch (sigma)
		case {'SM'}
			[void idx] = sort(abs(E1),1,'ascend');
			E1 = E1(idx);
			E = E1(1:k);

	end
end

end % eigs_lanczos


function [V E Q T] = remove_spurious(Q, T, sigma, k)
		if (idx >= k+ks_+1)
			[S Theta] = eigs(T(1:idx,1:idx), k+ks_, sigma);
			E_hat = sort(eigs(T(2:idx,2:idx),k+ks_,sigma));% min(idx-1,k+ks),sigma));
			[E p]     = sort(diag(Theta));

			q = 1; %q_mutiple = zeros(size(E));
			% remove multiple eigenvalues
			for jdx = 2:length(E)
				if (abs(E(jdx) - E(jdx-1)) < max(beta*(abs(S(idx,p(jdx-1)) + abs(S(idx,p(jdx))))),tol))
					% do not keep
					%q_mult(jdx) = 0;
				else
					% keep
					q = [q jdx];%$q_mult(jdx) = 1;
				end
			end
%			[[diff(E); 0] E-0.5 beta*S(end,:)']
			E = E(q);
			S = S(:,p(q));
%			[E-0.5 beta*S(end,:)']
			ks_new = max(0,k-length(E));

%{
			% remove spurious eigenvalues
			j_hat = 1;
			q_spur = [];
			for jdx=1:length(q)
				if (1 == q_mult(jdx))
					while (E_hat(jdx) < E(q(idx))
						j_hat = j_hat + 1;
					end
                        	        if (    ( j_hat > 1 && abs(E(jdx) - E(h_hat-1)) < 2*beta*(abs(S(idx,p(idx-1)))) ) ...
					     || ( j_hat < length(E)+1 && abs(E(jdx) - E_hat(j_hat)) < 2*beta*(abs(S(idx,p(idx-1)))) ) )
						% do not keep
					else
						q2 = [q2 q(jdx)];
					end
				end
			end		
% {
			% find spurious eigenvalues
			q = zeros(k+ks,1);
			ks = 0;
			l1 = 1;
			l2 = 1;
			lq = 0;
			ok = 1;
			while ( l1 <= length(E) && l2 <= length(E2) )
				if (abs(E(l1) - E2(l2)) < abs(E(l2))*sqrt(eps))
					% spurious eigenpair, do not keep
					ok = 0;
				end
				if (E(l1) < E2(l2))
					if (ok)
						lq = lq+1;
						q(lq) = l1;
					else
						ks = ks+1;
						ok =1;
					end
					l1 = l1+1;
				else
					l2 = l2+1;
				end			
			end
			q(lq+1:length(E)-l1) = l1+1:length(E);
			z = ones(length(E),1); z(q)=0;
			[[E2-0.5; inf] E-0.5 z]
%}
			%for jdx=2:length(E)
				%if (abs(E(jdx)-E(jdx-1)) < sqrt(eps))
				%	% spurious eigenpair, do not keep
				%	ks = ks + 1;
				%else
				%	% not spurious, keep eigenpair
				%	q = [q jdx];
				%end
			%end
			%ks = ks+2;
			%[E-0.5 abs(beta)*(S(end,:)./E')'.^2]
			%E = E(q(1:k));
			%S = S(:,p(q(1:k)));
	end
end % remove_spurious

function old_inner
%function [V E Q T R err_E err_Q] = eigs_lanczos(A, k, m, sigma)
	for idx=2:m-1
		idx
		k_ = min(k, idx-1);
		% compute k first eigenpairs of the projected sytems
		[S Theta] = eigs(T(1:idx,1:idx), k_, sigma);
		[E p]     = sort(diag(Theta));
		S = S(:,p);

		E2 = sort(eigs(T(2:idx,2:idx),k_, sigma));
		%[V E T Q] = remove_spurious(T, Q);
		%[V E T Q] = reorthogonalise(T, Q);
		
		% convergence estimate
		R(idx) = abs(beta(idx+1))*max(abs(S(end,:)./E'));

		% convergence statistic
		if (nargout > 5)
			%l = min(k-ks, length(E));
			l = k_;
			err_E(idx,1) = max(abs((E(1:l) - E_true(1:l))./E_true(1:l)));
			err_Q(idx,1) = norm(Q(:,1:idx)'*Q(:,1:idx)-speye(idx));
		end

		%if (length(E) >= k && R(idx) < tol)
		%	break;
		%end
		%if (idx >= maxiter)
		%	'error, No Convergence'
		%	break;
		%end
	end % for idx

	if (nargout < 2)
		% only return eigenvalues
		V = E(1:k,1);
	else
		% compute eigenvectors
		V = Q(:,1:end)*S(:,1:k);
		E = diag(E(1:k,1));
	end
end % old_inner

