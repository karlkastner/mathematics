% Wed Jul 11 18:52:29 MSK 2012
% Sun Jan  5 10:07:25 WIB 2014
% Karl Kästner, Berlin
%
% Thu Jan 22 11:13:23 CET 2015, nv changed to order = (nv+1)
% note, there was the parallel function vandermonde, for which n := n+1
%
%% van der Monde matrix
function A = vander_1d(x,order)
	if (~isscalar(order))
		x     = cvec(x);
		n     = sum(order~=0);
		order = find(order);
		A     = zeros(size(x,1),n);
		for idx=1:n
			A(:,idx) = x.^(order(idx)-1);
		end	
	else
		if (isvector(x))
			x = cvec(x);
		end
		A = vander_1d_(x,order).';
		s = size(x);
		A = reshape(A,[order+1,s(1),s(2)]);
		A = permute(A,[2, 1, 3]);
	end
end % vander_1d

function A = vander_1d_(x,order)
	n = size(x,1);
	A = zeros(n,order+1,class(x));
	A(:,1,:) = 1;
	if (order > 0)
		A(:,2) = x;
	end
	for idx=3:order+1
		A(:,idx) = A(:,idx-1).*x;
	end
end % vander_1d_

