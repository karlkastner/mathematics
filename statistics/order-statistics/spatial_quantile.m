% Sa 16. Jan 13:51:36 CET 2016
% Karl Kastner, Berlin
%% spatial quantile
function Q = spatial_quantile(X,p)
	m = size(X,1);
	R = spatial_rank(X);
	Q = zeros(m,length(p));
	p = [1-2*p];
	clf
	for idx=1:m
		[ri sdx] = sort(R(idx,:));
		plot(ri)
		hold on
		Q(idx,:) = interp1(ri,X(idx,sdx),p);
	end
	pause
end

