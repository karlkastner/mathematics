% Wed 15 May 14:31:34 CEST 2024
% Karl Kastner, Berlin
function q = skewinv(P,mu,sd,sk)
	q = zeros(size(P));
	for idx=1:numel(P)
		% initial guess: inverse of the normal distribution
		q0 = norminv(P(idx),mu,sd);
		q(idx) = fzero(@(q) skewcdf(q,mu,sd,sk) - P,q0);
	end % for idx
end

