% Tue Nov  1 00:50:04 MSK 2011

% derivation of end-point interpolation
% n = 1 : linear, f(x) = a1 x + a0
% n = 2 : quadratic, f(x) = a2 x^2 + a1 x + a0
% n = 4 : quartic, f(x) = a4 x^4 + a3 x^3 + a2 x^2 + a1 x + a0
% interpolation of last point is like interpolation of first point, but coefficients flipped

% m : outer point to interpolate - higher accuracies require more than one ghostpoints outside the domain
function derive_interpolation(n,m,vargrid)
	if (nargin < 2)
		m=1;
	end
	if (nargin < 3)
		vargrid = 0;
	end
	if (0 == vargrid)
		A=diag(0:n)*ones(n+1);
		A=A.^(A')
		s = (-m).^(0:n)
		s*inv(A)
		%symbolic:
		%for idx=0:n; F(idx+1,1)=sym(sprintf('f%d',idx)); end
		%H = sym('-h').^(0:n)
		%D = sym('h').^(0:n)
		%H*((A*diag(D)) \ F)
		%s*(A\F)
	else
		% variable coefficient
		P=ones(n+1)*diag(0:n)
		 %for idx=0:n; X(idx+1,1)=sym(sprintf('x%d',idx)); end
		 %A = X(1)-diag(X)*ones(n+1);
		 %A = A.^P
		 %H = (X(1)-X(2)).^(0:n)
		 %s = H *(A \ F)
		 %s = A \ F
		
		H=sym('0'); for idx=1:n; H(idx+1,1)=sym(sprintf('h%d',idx)); end
		H = diag(H)*ones(n+1);
		A = H.^P
		% choose hm1 to be h1
		% choose f0 to be zero
		H =  sym('-h1').^(0:n)
		H = [1 -h1 h1^2 -h1^3 h1^4]
		s=H*inv(A)
		s = simple(s)
	end % variable coefficients
end % derive_interpolation

