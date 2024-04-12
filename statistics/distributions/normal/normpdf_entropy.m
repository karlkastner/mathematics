% Tue 17 Jan 14:49:11 CET 2023
% Entropy of the normal distribution
% function h = normpdf_entropy(sd)
function h = normpdf_entropy(sd)
	if (issym(sd))
		pi_ = sym(pi);
	else
		pi_ = pi;
	end
	h  = 1/2*(1 + log2(2*pi*sd^2));
end

