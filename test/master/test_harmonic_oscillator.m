% Mon Feb 20 19:03:14 MSK 2012
% Karl Kästner, Berlin

function test_harmonic_oscillator()
%k = 11;
%N=2.^(4:10);
d=2
a=0.5*sqrt(2);
k = 11;
N=2.^(4:6);
E = zeros(k,length(N));
T = zeros(1,length(N));

for idx=1:length(N)
	tic()
	n=N(idx)
	A = harmonic_oscillator(n, d, a);
	E(:,idx) = eigs(A, k, 'SM')
	T(idx) = toc()
end

Err = E - E(:,end)*ones(1,size(E,2))

sqrt(sum(Err.^2))

end % function harmonic_osciallator()



