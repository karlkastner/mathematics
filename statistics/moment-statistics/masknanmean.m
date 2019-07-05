% 2019-06-27 15:05:29.778851052 +0200
%
% nanmean of the masked values of X
function m = masknanmean(X,mask)
	m = nanmean(X(mask));
end
