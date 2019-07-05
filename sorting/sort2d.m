% Mon Jun 24 20:03:34 UTC 2013
% Karl Kästner, Berlin
%
%% sort elements of matrix in X
%% returns row and column index of sorted values
%
function [v, id, jd] = sort2d(X)
	[v,i1] = sort(X(:));
	[id,jd]   = ind2sub(size(X),i1);
end
