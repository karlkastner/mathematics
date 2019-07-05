% Sun Jan  5 10:07:25 WIB 2014
% Karl Kastner, Berlin
%
%% van der monde matrix
function A = vandermonde(X, n)
	A = zeros(size(X,1),n+1);
	A(:,1) = 1;
	for idx=2:n+1
		A(:,idx) = A(:,idx-1).*X;
	end
end % vandermonde

