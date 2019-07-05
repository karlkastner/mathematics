% Fr 15. Jan 10:05:09 CET 2016
% Karl Kastner, Berlin
%
%% signed rank
%% Note: this is only a true rank if X is normal with zero mean, abitrary variance
function Q = spatial_singned_rank(X)
	Q = zeros(size(X));
	n = size(X,2);
	for idx=1:n
		%for jdx=[1:idx-1 idx+1:n]
		%	Q(:,idx) = Q(:,idx) + spatial_sign(X(:,idx)-X(:,jdx)) ...
		%		            + spatial_sign(X(:,idx)+X(:,jdx));
		%end
		Q(:,idx) = sum(spatial_sign(bsxfun(@minus,X(:,idx),X(:,[1:idx-1,idx+1:n]))),2) ...
		         + sum(spatial_sign(bsxfun(@plus,X(:,idx),X(:,[1:idx-1,idx+1:n]))),2);
	end
	Q = (0.5/(n-1))*Q;
	sum(Q,2)
	clf
	plot(sort(Q'))
	pause
end

