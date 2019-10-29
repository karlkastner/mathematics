% 0.5 * \Delta (1/r phi) + 1/r^2 phi = -lambda 1/r phi
% 0.5*(1/r \Delta phi + 2 \grad 1/r \grad phi + \Delta 1/r phi) + 1/r^2 phi = -lambda 1/r phi

L0 = 10;
odd = 0;
N=2.^(2:9) + odd;

k=10;
E_ = zeros(k,length(N));

for idx=1:length(N);
	n=N(idx)
	nh = ((n^2)+odd)/2
	h = 2*L0/(n+1);
	I = speye(n);
%=1;
	x = diag(sparse(2*L0*((1:n)/(n+1)-0.5)));
	X = kron(x,I); Y=kron(I,x);
	R = sqrt(X.^2 + Y.^2);
	Ri = diag(sparse(1)./diag(R));
if (odd)
	Ri(nh,nh) = 1e14; % quick fix
end
	D = 1/h*spdiags(ones(n,1)*0.5*[-1 0 1],-1:1,n,n);
	D = kron(I,D) + kron(D,I);
	L = 1/h^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	L = kron(I,L) + kron(I,L);
	D_Ri=D*Ri;
	%Riv = diag(Ri); D_Ri = diag(D*Riv);
	%D_Ri = diag(-(diag(X) + diag(Y)).*diag(Ri).^3);
	%full(diag(D_Ri_)) %(1:10,nh+1:nh+10))
	%full(diag(D_Ri)) %(nh+1:10,nh+1:10))
	%Derr = diag(diag(D_Ri)./diag(D_Ri_));
	%full(Derr(nh+1:nh+10,nh+1:nh+10))
	%Derr = full([2*diag(D_Ri_) diag(D_Ri) 2*diag(D_Ri_)-diag(D_Ri)])
	%	full(Derr)
	%norm(full(diag(D_Ri_)-diag(D_Ri)))/n^2
	a = 1e-1;
	V = diag(sparse(1./sqrt(diag(X).^2 + diag(Y).^2 + a^2*ones(n^2,1))));
	A = 0.5*(Ri*L + L*Ri) + D_Ri*D + V*Ri;
	% impose BC
if (odd)
	A = [A(1:nh-1, 1:nh-1) A(1:nh-1, nh+1:end)
	     A(nh+1:end, 1:nh-1) A(nh+1:end, nh+1:end)];
	Ri = [Ri(1:nh, 1:nh) Ri(1:nh, nh+2:end)
	     Ri(nh+2:end, 1:nh) Ri(nh+2:end, nh+2:end)];
end
%	norm(A-A')
	kk = min(size(A,1),k);
	E = eigs(A,-Ri,kk,'SM');
	E_(1:kk,idx) = sort(real(E))
end

E_
Err = E_ - E_(:,end)*ones(1,length(N))
sqrt(sum(Err.^2))
plot(E_(:,end))
E_ - E_old
E_old = E_;

