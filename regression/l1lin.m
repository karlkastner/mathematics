% 2014-08-25 19:25:18 +0200
% Karl Kastner, Berlin
%
%% solve ||Ax - b||_L1 by means of linear programming
%
function [x, E, X] = l1lin(A,b)
	ns = size(A,1);
	nc = size(A,2);
	%x  = randn(ns+nc,1);
	%f = @(x) sum(x(nc+(1:ns)))
	f  = [zeros(nc,1); zeros(ns,1)]';
	I = eye(ns);
%	Z = zeros(size(A));
%	AA = [A -I; -A -I; Z -I];
%	bb = [b;-b;zeros(size(b))];
%	AA = [A -I; -A -I];
%	bb = [b;-b];
	AA = [-A -I; A -I];
	bb = [-b;b];
	opt = optimset('TolFun',1e-14);
	lb = [];
	lb = [-inf*ones(nc,1); zeros(ns,1)];
	[x fval exitflag] = linprog(f,AA,bb,[],[],lb,[],[],opt);
%	fval
%	exitflag
%	f*x
	%x = linprog(f,AA,bb);
	s = x(nc+1:end);
	x = x(1:nc);
	err = A*x-b;
	E = norm(err);
	X = x;

s
	
	
	

if (0)
% from: Linear Programming: Foundations and Extensions
% TODO This algorithm does _not_ work
	reltol = 1e-7;
	maxiter = 100;
	p = 1;
	% take least squares solution as start point
	x = A \ b;
	x = randn(size(A,2),1);
	idx = 0;
	while (idx < maxiter)
		e = -(A*x - b);
		idx=idx+1;
		ne = sum(abs(e));
		E(idx) = ne;
		X(:,idx) = x;
		if (sum(abs(e)) < reltol*sum(abs(x)))
			return;
		end
		fdx = find(abs(e) > reltol*sum(abs(e)));
		fdx=(1:length(x))';
		ei = zeros(size(e));
		ei(fdx) = 1./e(fdx);
		%ei(fdx) = 1./(e(fdx)+sign(e(fdx))*ne*reltol);
		%x_ = (A'*diag(ei)*A) \ (A'*diag(ei)*b);
		x_ = (A(fdx,:)'*diag(ei(fdx))*A(fdx,:)) \ (A(fdx,:)'*diag(ei(fdx))*b(fdx));
		x = p*x_ + (1-p)*x;
		%dx = norm(x-x_);
		%x = x_;
		if (idx==maxiter)
			error
		end
	end
end
end

