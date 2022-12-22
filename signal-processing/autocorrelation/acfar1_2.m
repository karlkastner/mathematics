% 2017-08-11 22:02:46.822823392 +0200
%% autocorrelation of the ar1 process
function [a b] = acfar1_2(rho,n,id,flag)
	bfun = @(k)          [rho^k;
                            -2/(n*(n-k))*sum_ij(rho,n,n-k);
                            %-1/(n*(n-k))*sum_ij(rho,n,n-k);
		            % -1/(n*(n-k))*sum_ij_partial_(rho,k,n,n)
		            ... -1/(n*(n-k))*sum_ij_(rho,n,n-k);
		             mu2ar1(1,rho,n)]
	for idx=1:length(id)
		k = id(idx);
		b(idx,:) = (bfun(k));
	end
	b
	b = sum(b,2);
	if (nargin < 3 || nargin()>3 && flag)
		b = (n-id)./n.*b;
	end
	b0=sum(bfun(0));
	a = b/b0;
end

