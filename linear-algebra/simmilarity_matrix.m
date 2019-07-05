% 2013-08-06 23:46:57.000000000 +0200

function W = simmilarity_matrix(W,flag)
	W = distmat(W,flag);
	% distance to similarity
	fdx = find(W(:) > 0);
	W(fdx) = 1./W(fdx);
	W = W + speye(size(W));
end

