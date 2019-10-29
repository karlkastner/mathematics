% 2013 Aug 12  2012, Karl Kästner, Berlin
% Tue Mar 19 17:38:01 MSK 2013

%n=10; A = poisson([n n n]); E = eig(full(A)); A=A+1/12*A^2/(n+1)^2; H=eig(full(A)); [A B] = poisson_fem([n n n]); F = eig(full(A),full(B)); G = -pi^2*(length(F):-1:1)'.^2; plot([E F 0.5*(E+F) H]); %semilogy(abs([E-G 0.5*(E+F)-G F-G]))  
%n=100; A = poisson([n]); E = eig(full(A)); A=A+1/12*A^2/(n+1)^2; H=eig(full(A)); [A B] = poisson_fem([n]); F = eig(full(A),full(B)); G = -pi^2*(length(F):-1:1)'.^2; plot([E F 0.5*(E+F) G H]); %semilogy(abs([E-G 0.5*(E+F)-G F-G])) 
%n=100000; A = poisson(n); tic(); eigs(A,[],1,'SM'); toc(), [A B] = poisson_fem(n); tic(); eigs(A,B,1,'SM'); toc()
%n=100; A = poisson([n n]); tic(); eigs(A,[],10,'SM'); toc(), [A B] = poisson_fem([n n]); tic(); eigs(A,B,10,'SM'); toc()

function test_multigrid()
	opengl neverselect;
	N = 2.^(1:13) - 1; %unique(ceil(2.^(0:1:5)));
	
	for idx=1:length(N)
		n = N(idx)
		A = poisson([n]);  % 1D
		b = rand(size(A,1),1);
%		b = randn(size(A,1),1);
		x0 = zeros(size(A,1),1);
%		x0 = rand(size(A,1),1);
%		x_ = A \ b;
		tol = 1e-7;
	%       [x nrvec_1 r] = jacobi(A,b,x0,tol);
	%	M(idx,1) = length(nrvec_1);
		M(idx,1) = 0;
		[x nrvec_2] = mg(A,x0,b,tol);
		M(idx,2) = length(nrvec_2);
	end
	subplot(1,2,1)
%		semilogy(nrvec_1,'.-b'); hold on
		semilogy(nrvec_2,'.-r');
	subplot(1,2,2)
	semilogx(N,M,'.')
end


%function [x flag nrvec r] = v_cycle(Ac, Mc, b, x, tol)
%	
%	x = v_step(Ac, Mc, b_, x_, idx);
%		
%	end
%end

% todo, test for maxiter
function [x R] = mg(A,x,b,tol)
	k_max = length(x);
	k = 1;
%	y = A \ x;
	clf();
	while (1)
		r = A*x - b;
		R(k) = sqrt(r'*r);
		if (R(k) < tol)
			% convergence
			break;
		end
		k = k+1;
		if (k > k_max)
			'mg: no convergence'
			break;
		end
		opt.g = 0; % no galerkin
		opt.d = 1; % no damping
		opt.s = 0; % jacobi-method
%		x = step(A,x,b,d);
		x = x - step(A, r, @(n) contractor(n,1), opt);
%		subplot(2,1,1);
%		plot(x-y); hold on
	end
%	subplot(2,1,2);
end % mg

% Multigrid step to solve Ax = b 
% d    : damping
% f(n) : function to constract contractor / expansion matrix
% s    : smoothener 0 : Jacobi, 1 : Gauss-Seidel
function e = step(A, r, f_contractor,opt)
	n = length(r);
	if (1 == n)
		% terminate recursion
		e = r / A;
	else
		% A(x0 + e) = b + r
		% r  = A*x - b = Ae
		% x0 = 0 => r = -b

		% contractor/restrictor (expander/interpolator)
		[E C] = feval(f_contractor, n);
%		C = 0.5*E';

		if (0 == opt.g)
			% setup of the differential operator on the small grid
			% TODO, make a function call here
			m  = (n+1)/2-1;
			Ac = pois(m);
		else	% Galerkin approach
%			Ac = C*A*E;
%			Ac = C*A*C';
			Ac = E'*A*E;
		end
%			if(n<20)
%			full(Ac)
%			pause
%			end

		rc = C*r;

		% smooth contracted error
		ec = step(Ac, rc, f_contractor, opt);
		% apply smoothing of the contracted system
		e = E*ec;

		% smooth (single step)
		switch (opt.s)
			case { 0 }
				% jacobi
				% x = D^-1(b - Rx)
				e = opt.d*(r - tril(A,-1)*e - triu(A,+1)*e)./diag(A) + (1-opt.d)*e;
			case {99}
				% nothing (matlab bug if second case is omitted)
			default
				% lexicographic gauss-seidel
				% x = U^-1(b - Lx)
				e = opt.d*(triu(A) \ (r - tril(A,-1)*e)) + (1-opt.d)*e;		
		end % switch
	end % if not terminal step
end % function step

function [x R r] = jacobi(A,b,x,tol)
	tol = 1e-3;
	k_max = 10*length(x);
	k = 1;
	Di = 1./diag(A);
	while (1)
		r = A*x - b;
		R(k) = sqrt(r'*r);
		if (R(k) < tol)
			% convergence
			break;
		end
		k = k+1;
		if (k > k_max)
			'jacobi: no convergence'
			break;
		end
		x = Di.*(b - tril(A,-1)*x - triu(A,+1)*x);
	end
end % function jacobi

function [E C] = contractor(n,d)
	m = (n+1)/2;
	m_ = m-1;
	C = 1*sparse( 1:m_, 2:2:2*m_, ones(m_,1), m_, 2*m_+1);
	E = sparse( [2:1:m 1:1:m 2:1:m+1], [2:2:2*m-1 1:2:2*m-1 1:2:2*m-1], [ones(m-1,1); 0.5*ones(m,1); 0.5*ones(m,1)] );
	E = 0.5*E(2:end-1,:)'; % 0.5 to compensate that error is spread over sum([0.5 1 0.5])
	E = 1.7*E; % TODO why is this correction factor necessary ? why does it only work in combination with jacobi?

if (n < 10)
%	full(C)
%	full(E)
%	pause
end

%	C = sparse( [1:n/2 1:n/2], [2*(1:n/2)-1 2*(1:n/2)], ones(n,1) );
	switch (d)
		case { 1 }
			% nothing to do
		case { 2 }
			C = (1/2)*C;
			I = speye(n);
			C = kron(I,C) + kron(C,I);
		case { 3 }
			C = (1/3)*C;
			I = speye(n);
			C = kron(I,C) + kron(C,I);
		default 
			error('contractor')
	end % switch
end % function contractor

function Ac = pois(m)
	Ac = (m+1)^2*spdiags(ones(m,1)*[1 -2 1], -1:1, m, m);
end

