% 2015-08-05 17:07:19.829213427 +0200
% Karl Kastner, Berlin
%
%% Autocorrelation function of the finite AR1 process
%%
%% a_k = 1/(n-k)sum x_ix_i+1 + (xi + xi+k)mu + mu^2
%%     = r^k + 1/n sum_ij + 1/n
%
% function [a b] = acfar1(rho,n,last)
function [a b] = acfar1(rho,n,id,biased)
	if (nargin() < 3)
		% last = n-1;
		id = 0:min(100,n-1);
	end
	if (nargin() < 4)
		biased = true;
	end

	mu2 = mu2ar1(1,rho,n);

if (1)
	bfun =  @(k)    rho.^k ...
                         - 2./(n*(n-k)).*( (n-k).*(1+rho)/(1-rho) ...
                                          -rho/(1-rho)^2*(1-rho.^(n-k)+rho.^k - rho^n) ) ...
			+ mu2;
%	bfun =  @(k)    rho.^k ...
%                         - 1./(n*(n-k)).*( (n-k).*(1+rho)/(1-rho) ...
%                                          -2*rho/(1-rho)^2*(1-rho.^(n-k)+rho.^k - rho^n) ) ...
%                  	-2/n*1/(1-rho)^2*rho*(1-rho^n);

	b = bfun(id);
	b0 = bfun(0);
	if (biased)
		b = ((n-id)/n).*b;
	end
else
	bfun =  @(k) (n-k).*(mu2 + rho.^k) ...
                         - 1/n*( 2*(n-k).*(1+rho)/(1-rho) ...
                               - 2*rho/(1-rho)^2*(1-rho.^(n-k)+rho.^k-rho^n) );


%	b_ = b;
	b = bfun(id);
	if (~biased)
		b = (n./(n-id)).*b;
	end
	b0 = bfun(0);
%	norm(b-b_)
%%	pause
	
end
	
	a  = b/b0;
end

