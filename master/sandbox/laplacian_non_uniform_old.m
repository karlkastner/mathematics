% 2011 Dec  5 14:10 MSD
% Karl Kästner, Berlin

% function expects one more point outside the domain
function [D L] = laplacian_non_uniform(x)
	n  = length(x)-2;
	xl = x(1:end-2);
	xc = x(2:end-1);
	xr = x(3:end);

	K_num = [ (xr - xc), (xl - xr), (xc - xl) ];
	% todo - multiply out due to cancellation errors
	K_den = (xc - xr).*(xc - xl);
	L_den = (xl - xr);

	% symmetric laplacian
	L = 2*spdiags(K_num,-1:1,n,n);
	% diagonal scaling matrix
	D = diag(sparse(1./(L_den.*K_den)));

	if (nargout < 2)
		D = D*L;
	end
end % function laplacian_non_uniform

