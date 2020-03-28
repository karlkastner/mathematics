function test_laplacian_simple

%M = {'strong', 'lower', 'weak', 'none'}
M = {'strong', 'lower', 'central', 'none'}

l=100; n=1000;
for mdx=1:length(M)
mode = M{mdx}
k = 8;
%N = unique(ceil(logspace(1,3,100)));
N = 2.^(4:12)
%N = 10.^(1:4)
%order = 2;
order = 2*(1:5);
% analytic solution:
 %E_true = -pi^2*((k:-1:1)').^2; % laplacian
 E_true = 2*(1:k)'-1 % harmonic oscillator
for odx=1:length(order);
	%L = laplacian_simple(2*N(end),order(odx)); %order(end));
	%E_true = sort(eigs(L,k,'SM'))
	%E_(:,odx) = E_true;
	for ndx=1:length(N);
		n = N(ndx);
		L = laplacian_simple(n,order(odx),mode);
		V = (2*l*diag(sparse((1:n)'/(n+1)-0.5))).^2; % harmonic oscillator
		L = -L/(2*l)^2 + V;
		E = sort(eigs(L,k,'SM'));
		Err(ndx,odx) = norm(E_true - E);
		EE(:,ndx) = E-E_true;
	end
end
%Err
%E_
%E__ = 2*(1  - cos(2*pi*(1:k)'/N(end)))
subplot(2,2,mdx)
loglog(N,Err)
ylim([1e-12 1e3])
grid on
legend(num2str(order'))
title(mode)
end

end

