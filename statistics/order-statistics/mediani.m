% 2015-05-27 11:32:55.337733366 +0200
%
%% index of median, if median is not unique, any of the values is chosen
function [me, mdx] = mediani(X);
	if (isvector(X))
		X = cvec(X);
	end
	% TODO, interpolate when length of set is even
	m        = ceil(size(X,1)/2);
	[X, sdx] = sort(X);
	me   = X(m,:);
	mdx  = sdx(m,:);
end

