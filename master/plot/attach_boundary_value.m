function [X_ V_] = attach_boundary_value(X,V,L0,dimension)
	X_ = [-L0; X; L0];
	k = size(V,2);
	% attach boundary values (0)
	switch dimension
		case {1}
		V_ = [zeros(1,size(V,2)); V; zeros(1,size(V,2))];
		case {2}
			n = sqrt(size(V(:,1),1));
			for idx=1:k
				W = reshape(V(:,idx),n,n);
				W = [zeros(n+2,1) [ zeros(1,n); W; zeros(1,n) ] zeros(n+2,1)];
				W = reshape(W,(n+2)^2,1);
				V_(:,idx) = W;
			end
	end
end % function attach boundary value
