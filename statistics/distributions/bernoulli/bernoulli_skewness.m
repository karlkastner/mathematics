% 2024-05-17 12:45:12.891022319 +0200
% Karl Kastner, Berlin
function sk = bernoulli_skewness(p)
	sk = (1-2*p)./sqrt(p.*(1-p));
end

