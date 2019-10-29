% Fri Oct 28 12:09:19 MSK 2011
% Karl Kästner, Berlin
function x = mv(A, theta, u, x)
	% form product ( I - u*u')(A - theta I)(I - u*u')
	x = x - u*(u'*x);
	x = A*x - theta*x;
	x = x - u*(u'*x);
end

