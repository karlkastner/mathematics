% 2015-08-03 14:23:48.520248007 +0200
% Karl Kastner, Berlin

	clear;
	m = [1e4, 100];
%	a = linspace(0
	n = 100;
	L = 10;
	a = L*(1:n)'/n;
%cvec(linspace(0,L,n));
	for idx=1:n
		idx
		[x]       = skewrnd(a(idx),m(1),m(2));
		v = sqrt(mean(var(x)));
		x = x./v;
%		v = mean(var(x))
		sk(idx,1) = mean(skewness(x));
		sk(idx,2) = mean(qskew(x,0));
%		sk(idx,1) = idx/n;
	end
	sk = max(sk,0);
	plot(a,sk);
	%sk = [sk;-sk];
	%A=vander_1d(sk(:,2),3); % [0 1 0 1]),
	x = sk(:,2);
	x = log(x);
%	A = [x.^-1 x.^-0.5 x.^0 x.^0.5 x.^1 x.^1.5 x.^2];
	A = [x.^-1 x.^-0.5 x.^0 x.^0.5 x.^1 x.^1.5 x.^2 x.^2.5 x.^3];

	A\sk(:,1)

