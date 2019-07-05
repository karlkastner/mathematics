% 2015-07-15 20:51:43.068146538 +0200
% Karl Kastner, Berlin
%
%% generate random ar1 process
%% e1 = randar1(sigma,p,n,m)
%
function [x] = randar1(sigma,p,n,m,mu,dist)
	if (nargin() < 4)
		m = 1;
	end
	if (nargin() < 5)
		mu = 0;
	end
	if (nargin() < 6)
		dist = 'normal';
	end
	switch (dist)
	case {'normal'}
		e1 = sigma*randn(n,m);
	case {'gamma'}
		[a, b] = gamma_mode_to_parameter(mu,sigma);
		e1 = gamrnd(a,b,n,m)-mu;	
	case {'beta'}
		[a, b] = beta_mode_to_parameter(mu,sigma);
		e1 = betarnd(a,b,n,m)-mu;	
	otherwise
		error('randar1');
	end
	a = sqrt(1-p^2);
	s = a;

	% the variance of the first sample before filtering is scaled up,
	% otherwise is the variance of the first samples is underestimated
	% eq(1,:) = sigma \tilde eps = rho*0 + (1-rho)^(-1/2)*(1-rho)^2*sigma*\tilde eps
	e1(1,:) = e1(1,:)/a;

	% generate a correlated series
	% e1(idx,:) = a*e1(idx,:) + p*e1(idx-1,:);
	e1 = filter(a,[1 -p],e1);
%	for idx=2:n
%		e1(idx,:) = a*e1(idx,:) + p*e1(idx-1,:);
%	end
%	D = spdiags(ones(n,1)*[-p 1]/a,-1:0,n,n);
%	e1 = D\(e1);
	% apply mean
	x = e1 + mu;
end % randar1

