% 2016-07-13 09:19:45.892713951 +0200
%% norms of columns
function x = colnorm(x)
	x = sqrt(sum(x.^2));
end
