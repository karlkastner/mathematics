% 2024-05-17 13:30:58.231758570 +0200
% Karl Kastner, Berlin
function ku = bernoulli_kurtosis(p)
	ku = 3 + (1-6*(1-p).*p)./(p.*(1-p));
end
