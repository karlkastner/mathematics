% 2015-01-30 14:16:36.248945525 +0100
%
%% mean of the masked values of X
function m = maskmean(X,mask)
	m = mean(X(mask));
end

