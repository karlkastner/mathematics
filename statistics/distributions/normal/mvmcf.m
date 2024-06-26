% 2024-06-11 21:08:37.354554428 +0200
% Karl Kastner, Berlin
%
%% function c = mvncf(x,mu,S)
%% characteristic function of the normal distribution
function c = mvncf(x,mu,S)
	x = cvec(x);
	mu = cvec(mu)
	c = exp(1i*mu'*x + x'*S*x);
end

