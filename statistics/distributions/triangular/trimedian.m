% Thu 22 Mar 18:34:35 CET 2018
%
%% median of the triangular distribution
function m = trimedian(a,b,c)
	m   = a + sqrt((c-a)*(b-a)/2);
	fdx = b > (a+c)/2;
	m(fdx) = c - sqrt((c-a)*(c-b)/2);
end

