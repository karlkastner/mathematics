% 2015-03-04 15:28:32.732850327 +0100
% Karl Kastner, Berlin
%% invalidate values that exceed n-times the robust standard deviation
function x = filter_(x,n)
	q = quantile(x(:),[0.16, 0.5, 0.84]);
	m = q(2);
	s = 0.5*(q(3)-q(2));
	x( (x < m-n*s | x > m+n*s ) )= NaN(class(x));
end

