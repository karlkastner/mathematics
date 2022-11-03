% Sun 12 Jun 21:38:50 CEST 2022
function [B,l,D2] = brownian_noise_1d_laplacian(L,n,k)
	if (length(n)<2)
		n(2)=1;
	end
	if (nargin()<3)
		k = n;
	end
	k = min(k,n(1)-1);
	D2 = derivative_matrix_2_1d(n(1),L,2,'dirichlet');
%	TODO explicitly compute the eigenvectors, with the analytic expression
	if (1)
	[V,l]=eigs(D2,k,'sm');
	e   = randn(n(1),n(2));
	l   = -diag(l);
	lis = sqrt(1./l);
	if (~isfinite(lis(1)))
	lis(1) = 0;
	end
	B   = V*diag(lis)*V'*e;
	else
	B = sqrtm(full(-D2)) \ e;
	end
	B = B-B(1,:);
	% transform the brownian bridge into a brownian motion
	r = randn(1,n(2));
	x = sqrt(L)*(0:n(1)-1)'/n(1);
	B = B+x*r;
%full(D2)
%	(abs(full(V*diag(lis).^2*V'*D2)))
%pause
end
 
