% 15:05:40.873225890 +0200
%% extract 3x3 kernel of a 2D discretization matrix
function K = extract_kernel(A,n,id)
	% get element index in n(1)xn(2) matrix
	%id_ = (id(1)-1 + (id(2)-1)*n(2) ) + 1
	id_ = sub2ind(n,id(1),id(2));
	% row index
	rid_ = id_*ones(3);
	cid_ = id_ + [[-1,0,+1]-n(1)
		      -1,0,1;
		      [-1,0,+1]+n(1)];
	id_ = sub2ind(n(1)*n(2)*[1,1],rid_,cid_);
	K = A(id_(:));
	K = full(reshape(K,3,3));
end


