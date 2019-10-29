% 2013 Mar  8 02:58 Karl Kästner, Berlin
% 
% power method for complex or multiple first eigenvcector

% power method for symmetric
% power method for asymmetric -> slower convergence ?
% power method for duplicate eignvalues -> two different start values
% power method for complex eigenpairs -> 3 successive iterates are linear dependent, see below
% power method for smallest value -> shift by lambda max
% power method for second value -> find first and than deflate / orthogonalise to first vector

function [l1 l2] = powerc(A)
%A = diag([2 2 1]); %fails (sign), raleigh does work
A = diag([2 -2 1]); % fails (swaps complex with real part), raleigh does not work
A = [0 0 -2; 0 1 0; -2 0 0]; % succeeds, raileigh fails
	% unique largest : second 2x2 eigenvalue is zero, compare with raleigh quotient
	% two equal : second 2x2 eigenvalue is zero, compare with raleigh quotient
	% two equal with opposit sign, works

%A = [0 1; -2 1]
%A = [1 1 1; 0 0.1 0; -1 1 1]
n_ ='inf';
%n_=2;
x = rand(size(A,1),1)-0.5;
x=x/norm(x,n_);
k=100;
for idx=1:k
	x = A*x;
	norm(x)
	x = x/norm(x,n_);
	x
	pause
end
pq = [A*x x] \ (A*(A*x))	% rank deficient if not complex -> check this out first
%x1 = A*x;
%x2 = A*x1;
%n = norm(x2,n_);
%x = x/n;
%x1 = x1/n;
%x2 = x2/n;
%eig([x2 x1 x])
%pq = [x1 x] \ -x2
p = pq(1)
q = pq(2)
raleigh = x'*A*x/(x'*x)
l1 = - p/2 + sqrt(p^2/4 - q)
l2 = - p/2 - sqrt(p^2/4 - q)
x
end
