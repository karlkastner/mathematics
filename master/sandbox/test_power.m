function y = f(A,x)
	[V E] = eig(A);
	y = zeros(size(x));
	for idx=1:length(A)
		y = y + V(:,idx)*E(idx,idx)*V(:,idx)'*x;
	end
end
