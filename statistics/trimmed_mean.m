% Mi 13. Jan 14:16:05 CET 2016
% Karl Kastner, Berlin
%
%% trimmed mean
function mu = trmean(x)
	n = length(x);
	x = sort(x);
	if (n > 4)
		n4 = round(n/4);
		mu = mean(x(n4:n-n4));
	else
		mu = mean(x);
	end
end

