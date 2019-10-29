% Oct 22 23:56
% Karl Kästner, Berlin

% Modern methods in scientific computing and applications  By Anne Bourlioux, Martin J. Gander, Gert Sabidussi

function [Tv Tw D V W] = lanczos_biorthogonal_improved(A)
n=size(A,1);
% A*V = V*T_v
% A'*W = W*T_w
% W'AV = D*T_v, D = diag(d)
% V^-1 A V = T_v
%v = randn(n,1); %ones(n,1)/n;
%w = randn(n,1); %ones(n,1)/n;
v = zeros(n,1); v(1) = 1;
w = zeros(n,1); w(1) = 1;
rho = norm(v); v=v/rho;
nu = norm(w); w=w/nu;
d_old=1;
v_old = 0.*v;
w_old = 0.*w;
V(:,1) = v;
W(:,1) = w;
for idx=1:n
	%d = sqrt(w'*v); D(idx) = d;
	d = w'*v; D(idx) = d;
	a = w'*A*v/d;     AA(idx) = a;
	b = d/d_old*nu;   B(idx) = b;
	c = d/d_old*rho;  C(idx) = c;
	v_next = A*v - a*v - b*v_old;
	w_next = A'*w - a*w - c*w_old;
	Rho(idx) = norm(v);
	Nu(idx) = norm(w);
	d_old = d;
	v_old = v;
	w_old = w;
	v = v_next/Rho(idx);
	w = w_next/Nu(idx);


	if (idx < n)
		V(:,idx+1) = v;
		W(:,idx+1) = w;
	end
end % for idx

Tv = diag(B(2:end),1) + diag(AA) + diag(Rho(1:end-1),-1);
Tw = diag(C(2:end),1) + diag(AA) + diag(Nu(1:end-1),-1);

end % function l2

