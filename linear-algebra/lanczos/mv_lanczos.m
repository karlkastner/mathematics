% mv-lanczos
% computes y = T b - beta c
% with T = [s D s] (tridiagonal matrix with constant sub and superdiagonal)
% call: mv_1D(1/h^2,A(1:n+1:n^2),Q(:,idx), -beta, Q(:,idx-1))
function y = mv_1D(s,D,b)
	y = zeros(length(n));
	n = length(b);
	y(1) = D(1)*b(1) + s*b(2) + beta * c(idx)
	for idx=2:n-1
		y(idx) = s*b(idx-1) + D(idx)*b(idx) + s*b(idx+1) + beta * c(idx)
	end
	y(n) = s*b(n-1) + D(n)*b(n)
end

% two dimentsions
function 
	y = zeros(length(n));
	n = length(b);
	for jdx=2:m-1
		for idx=2:n-1
		end % for idx
		y(1) = s2 D(1)*b(1) + s*b(2) + beta * c(idx)
		for idx=2:n-1		
			y(idx) = s2*b(idx-n) + s1*b(idx-1) + D(idx)*b(idx) + s*b(idx+1) + s2*b(idx+n) + beta * c(idx)
		end % for idx
		y(n) = s*b(n-1) + D(n)*b(n)
		for idx=2:n-1
		end
	end % for jdx
end

