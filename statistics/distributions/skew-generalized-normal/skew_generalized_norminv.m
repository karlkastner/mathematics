% Wed 15 May 14:31:34 CEST 2024
% Karl Kastner, Berlin
function q = skew_generalized_norminv(P,mu,sd,l1,l2)
	q = zeros(size(P));
	for idx=1:numel(P)
		% initial guess: inverse of the normal distribution
		q0 = norminv(P(idx),mu,sd);
		q(idx) = fzero(@(q) skew_generalized_normcdf(q,mu,sd,l1,l2) - P,q0);
	end % for idx
end % skew_generalized_norminv 

