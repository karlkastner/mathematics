
%{
syms x1 x2 x3 y1 y2 y3 x y

%A = [ 1, x1, y1; 1, x2, y2;  1, x3, y3 ]
%C = inv(A).'*det(A) % det A for simplification
%A = [1 0 0; 1 1 0; 1 0 1]
A = [1 0 0; 1 1 0; 1 1 -1]
C = inv(A).'

phi = C*[1 x y].';

for jdx=1:3
for idx=jdx:3
	f = phi(idx)*phi(jdx);
	f = int(f,y);
	f = subs(f,y,(x));
	%f = subs(f,y,(1-x));
	f = int(f,x);
	f = subs(f,x,1);
	[n d] = rat(f);
	[idx jdx n d]
end
end
%}
%simple(f)
addpath('./int');
javaaddpath('.');
javaaddpath('/usr/share/java');
%[P T Bc] = mesh_2d_uniform([4 4],[3 3])
n=10;
% {
[P T Bc] = mesh_2d_uniform([n+2 n+2],[1 1]);
%P = P - 1;
mesh = Mesh(P,T,Bc);
mesh.prefetch_2d();

int = @int_2d_gauss_25; % 10th order
%int = @int_2d_trapezoidal;
%A = assemble_2d_dphi_dphi_java(mesh, [], int);
%B = assemble_2d_phi_phi_java(mesh, [], int);
A = assemble_2d_dphi_dphi(mesh, [], int);
B = assemble_2d_phi_phi(mesh, [], int);
f = @potential_coulomb;
V = assemble_2d_phi_phi(mesh, f, int);
[A B] = boundary_2d(A,B,Bc,1);
[v e] = eigs(A,B,1,'SM');
A = full(A);
B=full(B);
%[v e] = eig(A,B);
[e edx] = sort(diag(e));
v = v(:,end);
e = e(end);
% }
if (n<6)
	A
	12*B
end
	%full(V)

	x_ =  pi*[1:n]'/(n+1);
	I=ones(n,1);
	x=kron(I,x_);
	y=kron(x_,I);
	v_=sin(x).*sin(y);
v=-v;
%format long
e(end)
[A*v_ e(end)*B*v_]
%[A*v_ e(1)*B*v_]
v_'*A*v_/(v_'*B*v_)
v'*A*v/(v'*B*v)
v'*A*v/(v'*B*v)

	v_=v_/norm(v_);
%	d=3/(n+1)^3*(sin(x).*sin(y).*(cos(x).*cos(y)));
	d=3/(n+1)^3*(sin(x).*sin(y).*(cos(x).*cos(y).*(x-pi/2).^2.*(y-pi/2).^2));
	v=v(:,end);
	v=v/norm(v);
%[v_ v]
	v = v*sign(v(1))*sign(v_(1));
	norm(v-v_)
	subplot(2,3,1);
	plot([v-v_,d],'.-')
	subplot(2,3,2);
	imagesc(reshape(v,n,n));
	colorbar();
	subplot(2,3,3);
	imagesc(reshape(v_,n,n));
	colorbar();
	subplot(2,2,3);
	imagesc(reshape(v-v_,n,n));
	subplot(2,2,4);
	imagesc(reshape(d,n,n));

