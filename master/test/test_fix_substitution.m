% phi = r psi
% ( grad + 0.5 del + I) psi = lambda R psi
function test_fix_substitution()

N=2.^(4:16);
L = 100;
x0 = 0.5; % relative

for idx=1:length(N)
	n=N(idx);
	h = 2*L/(n+1);
	x = 2*L*((1:n)'/(n+1) - x0);
	R = diag(sparse(abs(x)));
	grad = 1/h*spdiags(ones(n,1)*0.5*[-1 0 1],-1:1,n,n);
	del  = 1/h^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	I    = speye(n);
	A1 = diag(sparse(sign(x)))*grad + 0.5*R*del + I;
	E__(:,idx) = sort(eigs(A1*R,R*R,6,'SM'));
	full(A1(1:10,1:10))
	A2 = 0.5*del*R + I;
	full(A2(1:10,1:10))
	%full(A1-A2)
	E(:,idx) = sort(eigs(A2,R,6,'SM'));
	snorm(A1-A2)
	A = 0.5*del + diag(sparse(1./abs(x)));
	full(A(1:10,1:10))
	A = lanczos(A, 100);
%	A
	E_(:,idx) = sort(eigs(A,[],6,'SM'));
%pause
end
E
D1 = E - E(:,end)*ones(1,length(N))
D1_ = sum(D1.*D1)
E_
D2=E_ - E_(:,end)*ones(1,length(N))
D2_ = sum(D2.*D2)
E__
D3=E__ - E__(:,end)*ones(1,length(N))
D3_ = sum(D3.*D3)
subplot(2,2,1)
loglog(N,[D1_;D2_;D3_]')
subplot(2,2,2)
loglog(N,abs(D1))
subplot(2,2,3)
loglog(N,abs(D2))
subplot(2,2,4)
loglog(N,abs(D3))
%E - E_
end

function n = snorm(A)
	n = sqrt(sum(sum(A.*A)));
end

