% 2015-08-05 11:49:36.846437699 +0200
% Karl Kastner, Berlin

s2 = [];
rho = 0.75;
N = (1:100)';
m = 1e4;
s = 2;
for idx=1:length(N)
	n = N(idx);
	X = randar1(1,rho,n,m);
	X = 2*X;
	%s2(idx,1) = varar1(rho,n);
	%s2(idx,2) = mean(var(X,1));
%	s2(idx,2) = mean(mean(X).^2);
	%s2(idx,2) = mean(var(X));
%	s2(idx,2) = mean(var(X));
%	s2(idx,3) = mean(mean(X).^2);

	s2(idx,1) = varar1(s,rho,n);
	s2(idx,2) = mean(var(X,1));


%	s2(idx,1) = mu2ar1(rho,n);
%	s2(idx,2) = var(mean(X));
end
s2
plot(N,s2)
