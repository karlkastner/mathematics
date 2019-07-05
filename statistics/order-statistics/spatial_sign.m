% 2016-01-15 09:35:15.496596195 +0100
% Karl Kastner, Berlin
%
%% spatial sign
function U = spatial_sign(X)
	% U = 1/||x||*X
	no = sqrt(sum(X.*X));
	U = bsxfun(@times,1./no,X);
end

