% Fri Dec 28 22:41:47 MSK 2012
% Karl Kästner, Berlin

% see HS 3.1a-3.1f

function [x nerrs a b] = cg(A,b)
	n = length(b);
%	x = zeros(n,1);
	x = rand(n,1);
	r = b - A*x;
	p = r;
	a = zeros(1,n);
	b = zeros(1,n);
	nerrs = zeros(1,n+1);
	nerrs(1) = r'*r;
	for idx=1:n
		Ap = A*p;
		a(idx) = nerrs(idx) / (p'*Ap);
		x = x + a(idx)*p;
		r = r - a(idx)*Ap;
		nerrs(idx+1) = r'*r;
		b(idx) = nerrs(idx+1)/nerrs(idx);
		p = r + b(idx)*p;
	end % for idx
end % function cg

