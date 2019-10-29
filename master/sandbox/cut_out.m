% 2012 Dec 19 16:59 MSK
% Karl Kästner, Berlin

% TODO: test for harmonic oscillator

% set-up of a 2D-laplacian matrix with cut-out region
r0 = 0.0;
N=2.^(4:8);
%n=100;

for idx=1:length(N)
	n = N(idx)
	A = laplacian_cut_out(n,r0);
	%A = A + speye(size(A));
	condest(A)
	E_ = eigs(A, 10, 'SM');
	[void id] = sort(real(E_));
	E(:,idx) = E_(id(:)) 
end

E
Err = E - E(:,end)*ones(1,length(N))

sqrt(sum(abs(Err).^2))
