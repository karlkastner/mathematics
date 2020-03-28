% Mon Oct 31 23:41:24 MSK 2011
% Karl Kästner, Berlin

% shifting indices makes no difference
% use format rat

sign inverted!! call with 2,0

function [dF M] = derive_laplacian(n, vargrid)
	if (nargin < 2)
		vargrid = 0;
	end
	if (0 == vargrid)
		C = [-n:-1 1:n]'*ones(1,2*n);
		D = diag(1./factorial(1:2*n))
		P = ones(2*n)*diag(1:2*n)
		A = (C.^P)*D
		B=[eye(n) -ones(n,1) zeros(n);
                    zeros(n) -ones(n,1) eye(n)]
		dF = A\B
		M = 1./(abs(dF(:,1)))
		dF = dF .* repmat(M,1,size(dF,2))
	else
	
	% with variable stepwidth
	for n=2
	% declare symbolic variables
	c=[-n:n];
	for idx=1:2*n+1
	if (c(idx) < 0)
		vstr{idx} = sprintf('xm%d',abs(c(idx)));
	elseif (c(idx) > 0)
		vstr{idx} = sprintf('xp%d',abs(c(idx)));
	else
		vstr{idx} = 'x0';
	end
	end
	for idx=1:n
	for jdx=1:2*n
		s = sprintf('((%s - %s)^%d)/%d', vstr{n+1}, vstr{idx}, jdx, factorial(jdx));
		A(idx,jdx) = sym(s);
		s = sprintf('((%s - %s)^%d)/%d', vstr{n+1}, vstr{n+1+idx}, jdx, factorial(jdx));
		A(n+idx,jdx) = sym(s);
	end
	end
	B=zeros(2*n,2*n+1);
	B(:,n+1) = -1;
	for idx=1:n
	B(idx,idx) = 1;
	B(2*n+1-idx,2*n+2-idx) = 1;
	end
	A
	B
	% A H dF = B F
	dF = A \ B;
	 %H = diag([h h^2 h^3 h^4]);
	 %F = [fm2 fm1 f0 fp1 fp2].';
	 %(A*H) \ (B*F)
	[num den] = numden(dF);
	num_ = factor(num)
	den_ = factor(den)
end

end % if vargrid

end % function derive laplacian

