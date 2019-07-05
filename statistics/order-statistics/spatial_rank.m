% Sa 16. Jan 13:50:03 CET 2016
% Karl Kastner, Berlin
%
%% unsigned rank
%
function R = spatial_rank(X)
	n = size(X,2);
	m = size(X,1);
	R = zeros(m,n);
	for idx=1:n
		R(:,idx) = sum(spatial_sign(bsxfun(@minus,X(:,idx),X(:,[1:idx-1,idx+1:n]))),2);
	end
	R = R/(n-1);
	sum(R,2)
end
