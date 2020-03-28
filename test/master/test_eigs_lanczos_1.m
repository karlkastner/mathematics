% Jul 13 21:49
% Karl Kästner, Berlin

function tel()

clear; figure(1); clf
L=10.0;
%N=[100 300 1000 3000];
%N = [200 2000];
N = [1000];
k = 10;
mode = 'SA';	% SA

for ndx=1:2 %length(N)
	n=N(ndx);
%	A = -0.5/(L/n)^2*spdiags(ones(n,1)*[1 -2 1], -1:1, n,n)-diag(sparse(n./(L*(1:n))));
%	A = -(n+1)^2*spdiags(ones(n,1)*[1 -2 1], -1:1, n, n); % + diag(sparse(1:n)).^2;
	%A = spdiags(ones(n,1)*[1 -2 1], -1:1, n, n); % + diag(sparse(1:n)).^2;
	n = 20;
	A1 = (n+1)^2*spdiags(ones(n,1)*[1 -2 1], -1:1, n, n); % + diag(sparse(1:n)).^2;
	I1 = speye(n);
	n=50
	A2 = (n+1)^2*spdiags(ones(n,1)*[1 -2 1], -1:1, n, n); % + diag(sparse(1:n)).^2;
	I2 = speye(n);
	n = 1000;
	A = kron(A1, I2) + kron(I1, A2);

	if (1) % == ndx)
%		[V_ E_] = eigs(A,k,mode);
%		x = sum(V_,2);
		x = rand(n,1); % - 0.5;
%		plot(x)
%		pause
	else
		x = interp1((0:N(ndx-1)-1)'/(N(ndx-1)-1), sum(V,2),(0:n-1)'/(n-1));
	end
	%[V E R err_E err_Q beta] = eigs_lanczos(A, k, mode, x, eigs(A, k, 'SM'));
	%[V E A_] = eigs_preconditioned(@eigs_lanczos, A, -0.5, x, k, mode);
	E_true = sort(eigs(A, k, mode));
	[V E Q T R err_E err_Q] = eigs_lanczos(A, k, mode);
	%[diag(E) - 0.5 sort(eigs(A, k, 'SM'))]
	[diag(E), E_true]

	subplot(2,2,ndx)
	%semilogy([R err_E err_Q beta(1:length(R))])
	semilogy([R err_E err_Q ],'LineWidth',2)
	grid on
	title(sprintf('Domain Radius %d a_0, Matrix Rank n = %d, Number of computed eigenvalues k = %d', L, n, k));
	legend('est. max |\lambda - \theta|', 'max |\lambda - \theta|', '|| Q''Q - I ||_2');
	xlabel('iteration i')
	ylabel('[1]')
	pause(1)
	ylim([1e-7 1e1]) 
end
print -depsc ../img/test_eigs_lanczos.eps

end

function [V E A] = eigs_preconditioned(solver, A, sigma, x, k, mode)
	A = A - sigma*speye(size(A));
	L = chol(A);
	A = inv(L)*A*inv(L)';
	% restore symmetrie
	A = 0.5*(A + A');
	[V E] = feval(solver, @(x) Afunc(A, L, x), sigma, k, mode)
end % eigs_preconditioned

function x = Afunc(Aml, L, x)
	x = L \ (Aml * ( L' \ x));
	%x = Aml*x;
end

