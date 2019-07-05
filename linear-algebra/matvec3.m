% Thu 11 Jan 14:40:19 CET 2018
% length of third dimension of A and b must either be 1 or identical
%
%% matrix-vector product of stacked matrices and vectors
function c = matvec3(A,b)
	% TODO allow for different matices and different b
	n = size(A,1);
	m = size(A,2);

	c = zeros(1,n,max(size(A,3),size(b,2)));

	b = shiftdim(b,-1);
	
	for idx=1:n
	 for jdx=1:m
		c(1,idx,:) = c(1,idx,:) + A(idx,jdx,:).*b(1,jdx,:);
	 end % for jdx
	end % for idx
	size(c)
%	c = shiftdim(c,-1);
	% was chnaged from -1 to plus 1 for tri_excircle
	c = shiftdim(c,+1);
end % matvec

