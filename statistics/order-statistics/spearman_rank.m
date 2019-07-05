% So 2. Aug 11:53:54 CEST 2015
% Karl Kastner, Berlin
%
%
function rho = spearman_rank(rx,ry,tx,ty)
	n = length(rx);
	rho = (n^3 - n - 0.5*tx - 0.5*ty - 6*sum((rx-ry).^2)) ./ ...
		(sqrt(n^3 -n -tx)*(n^3-n-ty));
	%rho = 1 - 6*sum((rx-ry).^2)/(n*(n^2-1));
end % spearman_rank

