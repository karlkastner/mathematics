% Fri Oct 28 12:09:19 MSK 2011
% Karl Kästner, Berlin

% stepwise form the product ( I - u*u')(A - theta I)(I - u*u')x
function x = mv_jacobi_davidson(A, theta, u, x)
	x = x - u*(u'*x);
	x = A*x - theta*x;
	x = x - u*(u'*x);
end % function mv_jacobi_davidson

