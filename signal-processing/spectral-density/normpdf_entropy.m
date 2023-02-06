% Tue 17 Jan 14:49:11 CET 2023
% Entropy of the normal distribution
function h = normpdf_entropy(sd)
	h  = 1/2*(1 + log2(2*pi*sd^2));
end

