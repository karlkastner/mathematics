function A = bernstein_matrix(x,n)
	nx = length(x);
	A = zeros(nx,n+1);
	for idx=1:n+1
		A(:,idx) = nchoosek(n,idx-1)*x.^(idx-1).*(1-x).^(n-idx+1);
	end
	
end

