% Fri  3 Jan 15:17:20 +08 2020
%% weighted geometric mean
%% function mu = wgeomean(w,x)
function mu = wgeomean(w,x)
	% mu = (prod x).^(1/n)
	%    = exp(1/n*sum(log(x))
	w = w/sum(w);
	mu = exp(sum(w.*log(x)));
end


