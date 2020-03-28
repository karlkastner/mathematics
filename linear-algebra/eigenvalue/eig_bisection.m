% Nov  9 19:39 MSC
% Karl Kästner

% n=4; L=diag(ones(n-1,1),-1) + -2*diag(ones(n,1)) + diag(ones(n-1,1),+1),
% I=eye(n); L2=kron(kron(L,I),I) + kron(kron(I,L),I) + kron(I,kron(I,L));
% H=hess(L2); D=diag(H,-1); Hs=hess(single(L2)); Ds=diag(Hs,-1); [abs(D)<1e-7 D
% abs(Ds)<1e-3 Ds]; sum(abs(D)<sqrt(eps))

% todo, split first in sections

function E = eig_bisection(D0,D1,k)
	n = length(D0);
	k = min(k,n);
	E = zeros(k,1);

	% lower bound on smallest eigenvalue
	e_l = min(D0 - [0; D1(2:end)] - [D1(2:end); 0]);
	% upper bound on largest eigenvalue
	e_u = max(D0 + [0; D1(2:end)] + [D1(2:end); 0]);

	for idx=1:k
		e = eig_bisection__(D0, D1, e_l, e_u, idx);
		E(k-idx+1,1) = e;
	end
end % eig_bisection

function e = eig_bisection__(D0, D1, e_l, e_u, k)
		tol = eps;
		n = length(D0);
		% evaluate polynomial
		while ( abs(e_l - e_u) > tol*(abs(e_l) + abs(e_u)) && abs(e_l) + abs(e_u) > tol)
			% central-point
			e_c = 0.5*(e_l + e_u);
			[p_c n_s] = tridet(D0-e_c, D1);
			% select half with sign change
			if (n_s > n-k)
				% sign change in lower half
				e_u = e_c;
			else
				% no sign change in lower half,
				% viz. sign change in upper half
				e_l = e_c;
			end
		end
		e = 0.5*(e_l + e_u);
end % eig_bisection_

function e = eig_bisection_(D0, D1, e_l, e_u, idx)
		tol = eps;
		% evaluate polynomial
		p_l = tridet(D0-e_l, D1);
		while ( abs(e_l - e_u) > tol*(abs(e_l) + abs(e_u)) )
			% central-point
			e_c = 0.5*(e_l + e_u);
			p_c = tridet(D0-e_c, D1);
			% select half with sign change
			if (p_l*p_c < 0)
				% sign change in lower half
				e_u = e_c;
			else
				% no sign change in lower half,
				% viz. sign change in upper half
				e_l = e_c;
				p_l = p_c;
			end
		end
		e = 0.5*(e_l + e_u);
end % eig_bisection_

% p   : value of determinant
% n_s : number of sign changes in sub-determinants
function [p n_s] = tridet(D0, D1)
	n = length(D0);
	n_s = 0;
	p_2 = 0;
	p_1 = 1;
	for idx=1:n
		p = D0(idx)*p_1 - D1(idx)^2*p_2;
		n_s = n_s + (p*p_1 < 0 || 0 == p);
		p_2 = p_1;
		p_1 = p;
	end
end % function pol

function [D0, D1, n] = deflate(e, D0, D1, n)
	% deflate the matrix
	A = spdiags([D1 D0 D1], -1:1, n,n);
	z = zeros(n,1);
	I = speye(n);
	% find eigenvector
	% todo trisolve
	v = (A - e*I) \ z;
	v
	full((A*v)./v)
	pause
	
	V = [v I(:,2:end)];
	A = V*A*V';
	full(A)
	pause
%	pause
%{
	D0 = diag(A);
	D0 = D0(2:end);
	size(D0)
	D0 = diag(A,+1);
	D1 = D1(2:end);
	size(D1)
	n = n-1;
%}
end % function deflate

