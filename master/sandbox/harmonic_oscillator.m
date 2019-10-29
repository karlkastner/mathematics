% Mon Feb 20 19:47:59 MSK 2012
% Karl Kästner, Berlin

%N=2.^(2:10); for idx=1:length(N); n=N(idx); d=1; A =  harmonic_oscillator(n, d); E(1:4,idx) = sort(eigs(A,[],4,'SM')), end


function [A X_] = harmonic_oscillator(n, d, a, b)
	L0 = 5;
	h = 2*L0/(n+1);

	X_ = h*(1:n)' - L0;
	X = diag(sparse(X_));
	I = speye(n);

	L  = 1/h^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);

	switch (d)
		case {1}
			A = (-L + X.^2);
		case {2}
			Y = X;
			A = -kron(I,L) + -kron(L,I) + kron(X.^2,I) + kron(I,a*Y.^2);
		case {3}
			Y = a*X;
			Z = b*X;
			A =     (  -kron(kron(L,I),I) ...
				 + -kron(kron(I,L),I) ...
				 + -kron(kron(I,I),L) ...
				 +  kron(kron(X,I),I).^2 ...
				 +  kron(kron(I,Y),I).^2 ...
				 +  kron(kron(I,I),Z).^2 );
		otherwise
			error
	end	
end

