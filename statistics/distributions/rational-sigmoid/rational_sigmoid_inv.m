% 2024-05-17 09:04:01.078334971 +0200
% Karl Kastner, Berlin
function x = rational_sigmoid_inv(F,mu,a,k)
	x = zeros(size(F));
	for idx=1:numel(F)
		x(idx) = fzero(@(x) F - rational_sigmoid_cdf(x,mu,a,k),0);
	end
end

