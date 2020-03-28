% Thu Mar  1 23:19:52 MSK 2012
% Karl Kästner, Berlin

function confinement()

d = 3;
eig_mode = 'SM'
k = 20;
ke = 20;
n = 2.^7;

H = 4.360e-18;
T = 1e3;
kb = 1.381e-23;


f = @hydrogen_2d

L = fliplr(linspace(1,2,11))%fliplr([1 2 4 8 16 32]);
eig_mode = 'SM' %SA
s = -2;
e=1/4; %^2

E = zeros(k,length(L));

for idx=1:length(L)
	L0 = L(idx);
	[A I] = feval(f,n,L0,e);
	A = 0.5*(A+A');
	E(1:min(k,n^d),idx) = sort(real(eigs(A-s*I,[],min(k,n^d),eig_mode)+s));
	E(1:ke,:)
end % for idx

log(sum(exp(-E*H/(T*kb))))

end % function fdm_vargrid

% todo, use setup vargrid, not laplacian non unigform
function [A I X L] = hydrogen_1d(n, L0,e)
	I = speye(n);
	h = L0/(n+1);
	X = h*(0:n+1) - L0/2;
	p = 4;	% 3 for 1d, exp for 2d, 4 for 3d
	c=1;
%	X = 0.5*L0/X(end)^p*X.^p.*sign(X.^(p-1));
	d = X(end)/(e*(exp(e*X(end))-1)); X = d*e*sign(X).*(exp(e*abs(X)) - 1);
	[L D1] = laplacian_non_uniform_2(X');
	L = sqrt(D1)*L*sqrt(D1);
	X = X(2:end-1);
	R = diag(sparse(abs(X)));
%	[L] = laplacian_non_uniform_2(X');
	A = -(0.5*L + inv(R));
end % hydrogen_1d

function [A2 I2] = hydrogen_2d(n,L0,e)
	[A1 I1 X L1] = hydrogen_1d(n, L0,e);
	L2 = kron(L1,I1) + kron(I1,L1);
	R2 = sqrt(kron(diag(sparse(X.^2)),I1) + kron(I1,diag(sparse(X.^2))));
	A2 = -(0.5*L2 + inv(R2));
	I2= speye(n^2);
end

function [A3 I3] = hydrogen_3d(n,L0,e)
	[A1 I1 X L] = hydrogen_1d(n, L0,e);
	L3 = kron(L,kron(I1,I1)) + kron(kron(I1,L),I1) + kron(I1,kron(I1,L));
	R3 = sqrt(kron(kron(diag(sparse(X.^2)),I1),I1) + kron(kron(I1,diag(sparse(X.^2))),I1) + kron(I1,kron(I1,diag(sparse(X.^2)))));
	A3 = -(0.5*L3 + inv(R3));
	I3= speye(n^3);
end

