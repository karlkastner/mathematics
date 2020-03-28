% radial schrödinger equation

function test_radial_fixes()

L0=100;
%N=1000;
%N=2.^(4:12);
N=2^12;
k = 10;
p = 2;

A_ = [0];
%A_ = [0 0.25 0.5 1]
%A_ = 1e-3;
A_ = [0 1e-4 1e-3 1e-2 1e-1 1]

E_true = -1./(1:10)';

for adx=1:length(A_)
	a = A_(adx);
	E = [];
for idx=1:length(N)
	n = N(idx);
	h = L0/(n+1);
	x=L0*(1:n)/(n+1);
	% potential
	Ri = diag(sparse(1./x));
	V = diag(sparse(1./(x.^p + a^p).^(1/p)));
	% laplacian
	L = spdiags(ones(n,1)*[1 -2 1],-1:1, n, n);
	%L(1,1:5) = L(1,1:5) + [5.0000  -10.0000   10.0000   -5.0000    1.0000];
	%L(1,1:4) = L(1,1:4) + [ 4.0000   -6.0000    4.0000   -1.0000];
	%L(1,1:3) = L(1,1:3) + [3 -3 1];
	%L(1,1:2) = L(1,1:2) + [2 -1];
	L = 1/h^2*L;
	% radial schrödinger
	%A = -(0.5*L + V);
	A = -(L + 2*V);
	% solve equation
	% eigs returns 0 if SA choosen and number of wanted eigenvalues is to small
	%E = sort(eigs(A,20,'SA'));
	% => choose SM and and pos_max > |pos_min|
	E = sort(eigs(A,50,'SM'));
end
	EE(1:10,adx) = E(1:10)
end
EE
Err = EE - E_true*ones(1,size(EE,2))
Err_ = EE - EE(:,1)*ones(1,size(EE,2));
A_
sqrt(sum(Err.^2))
sqrt(sum(Err_.^2))

end

