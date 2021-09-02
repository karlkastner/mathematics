% Sun 11 Jul 12:19:45 CEST 2021
% TODO allow for dirichlet and neumann boundary conditions
function M = kernel2matrix(n,K);
	nn  = prod(n);
	% allowcate memory
	buf = zeros(numel(K)*nn,3);
	% row index
	r = repmat((1:n(1))',n(2),1);
	% column index
	c = flat(repmat((1:n(2)),n(1),1));
	k = 0;
	sk = size(K);
	if (any(mod(sk,2)~=1))
		error('kernel must be odd');
	end
	sk = (sk-1)/2;
	for di=-sk(1):sk(1)
	 for dj=-sk(2):sk(2)
		% circular bc
		r_ = r+di;
		fdx = r_<1;
		r_(fdx) = r_(fdx) + n(1);
		fdx = r_>n(1);
		r_(fdx) = r_(fdx)-n(1);
		c_ = c+dj;
		fdx = c_<1;
		c_(fdx) = c_(fdx)+n(2);
		fdx = (c_>n(2));
		c_(c_>n(2))=c_(fdx)-n(2);
		buf(r+(c-1)*n(1)+k,1:3) = [r+(c-1)*n(1),r_+(c_-1)*n(1), K(di+sk(1)+1,dj+sk(2)+1)*ones(nn,1)];
		k = k+nn;
	 end
	end
	M = sparse(buf(:,1),buf(:,2),buf(:,3),nn,nn);
end

