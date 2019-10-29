% 2012 Jul 29 21:08
% Karl Kästner, Berlin

% 2D
function test_uniform_mesh(printflag)

format compact

N_ = 2.^(2:9)';
k=1;
L=20;
s = -2;
for idx=1:length(N_)
	n = N_(idx,1)

	% grid
	h = L/(n+1);
	%X = mesh_1d_uniform(n+2,L) - L/2;
	X = L*(-0.5:n+0.5)'/n - L/2;

	% analytic solution
	x = 2*X(2:end-1)*ones(1,n);
	y = x';
	w = 0.5*pi*exp(-sqrt(y.^2+x.^2));
	w4x=(-3.*y.^2.*(2.*x.^4 + 2.*x.^2.*y.^2 + 4.*x.^2 - y.^2)./(x.^2+y.^2).^(7/2) + (x.^6 + x.^4.*y.^2 - 12.*x.^2.*y.^2 + 3.*y.^4)./(x.^2+y.^2).^3).*w;
	w4y=(-3.*x.^2.*(2.*y.^4 + 2.*y.^2.*x.^2 + 4.*y.^2 - x.^2)./(y.^2+x.^2).^(7/2) + (y.^6 + y.^4.*x.^2 - 12.*y.^2.*x.^2 + 3.*x.^4)./(y.^2+x.^2).^3).*w;
	w4 = w4x+w4y;
	v_true  = reshape(w,n^2,1);
	v4_true = reshape(w4,n^2,1);
	e_true = -2;

	% FDM discretisation
	% difference operators
	D1 = 1/h*spdiags(ones(n,1)*[-1 0 1],-1:1,n,n);
	D2 = 1/h^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	D2_ = D2 - 1/12*h^2*D2^2 + 1/90*h^4*D2^3 - 1/560*h^6*D2^4;
	I = speye(n);
	D1 = kron(I,D1) + kron(D1,I);
	D2 = kron(I,D2) + kron(D2,I);
	D2_ = kron(I,D2_) + kron(D2_,I);
	D4 = D2^2;
	D4_ = D2_^2;
	I = speye(n^2);
	% potential
	XX = kron(ones(n,1),X(2:end-1));
	YY = kron(X(2:end-1),ones(n,1));
	V = diag(sparse(1./sqrt(XX.^2+YY.^2)));
	% Schrödinger PDE
	A = -0.5*D2 - V;
	% solve the eigenvalue problem
	[v_fdm e_fdm] = eigs(A - s*I,[],1,'SM');
	e_fdm = e_fdm + s;
	N(idx,1) = size(A,1);

	% FEM with linear polynomial basis functions
	% generate mesh
	[P T Bc] = mesh_2d_uniform([n n],[L L]);
	P = P-L/2;
	mesh = Mesh_2d(P,T,Bc);
	mesh.prefetch();
	% assemble the discretisation matrices
	f = Potential_2D_Coulomb;
	int = 'int_2d_gauss_3';
	A = assemble_2d_dphi_dphi_java(mesh,[],int);
	B = assemble_2d_phi_phi_java(mesh,[],int);
	V = assemble_2d_phi_phi_java(mesh,f,int);
	% Schrödinger PDE
	A = -0.5*A+V;
	% apply boundary conditions
	[A B p] = boundary_2d(A,B,Bc,1);
	% solve eigenvalue problem
	[v_fem_1_ e_fem_1] = eigs(A-s*B,B,1,'SM');
	% padd boundary values
	v_fem_1 = zeros(size(P,1),1); v_fem_1(p,1) = v_fem_1_;
	% undo shift
	e_fem_1 = e_fem_1 + s;
	N(idx,2) = size(A,1);

	% FEM with quadratic polynomials
	if (n<512)
	int = 'int_2d_gauss_6';
	mesh.promote_3_6();
	mesh.prefetch();
	A = assemble_2d_dphi_dphi_java(mesh,[],int);
	B = assemble_2d_phi_phi_java(mesh,[],int);
	V = assemble_2d_phi_phi_java(mesh,f,int);
	% Schrödinger PDE
	A = -0.5*A+V;
	% apply boundary conditions
	[A B p] = boundary_2d(A,B,Bc,1);
	% solve eigenvalue problem
	[v_fem_2 e_fem_2] = eigs(A-s*B,B,k,'SM');
	e_fem_2 = e_fem_2 + s;
	N(idx,3) = size(A,1);
	else
		e_fem_2 = NaN;
		v_fem_2 = NaN;
	end

	E(idx,:) = [e_fdm e_fem_1 e_fem_2];
	% eigenvalue error
	Err_e(idx,1)  = (e_fdm - e_true )/e_true;
	Err_e(idx,2)  = (e_fem_1 - e_true )/e_true;
	Err_e(idx,3)  = (e_fem_2 - e_true )/e_true;
	leg_e = {'FDM', 'FEM p=1, Gauss', 'FEM p=2, Gauss'}
	% eigenvector error
	Err_v(idx,1)  = norm(abs(v_fdm/norm(v_fdm)) - abs(v_true/norm(v_true)));
	Err_v(idx,2)  = norm(abs(v_fem_1/norm(v_fem_1)) - abs(v_true/norm(v_true)));
	Err_v(idx,3) = NaN;
	%Err_v(idx,2)  = norm(v_fem_2 - v_true);
	leg_v = {'FDM', 'FEM p=1, Gauss', 'FEM p=2, Gauss'}
	% error estimates
	Err_est(idx,1)  = h^2/(12*sqrt(v_true.'*v_true))*norm(v4_true,'inf');
	Err_est(idx,2)  = h^2/(12*v_fdm.'*v_true)*v_fdm.'*v4_true;
%	Err_est(idx,3)  = h^2/(12*v_true.'*v_true)*norm(D4*v_true,2);
	leg_est = {{'$\frac{h^2}{12 \sqrt{u_*^T u_*}} ||D_*^4 u_*||_\infty$', ...
	           '$\frac{h^2}{12 u_h^T u_*} |u_h D_*^4 u_*|$'}, 'interpreter', 'latex' }

%{		
	% finite difference derivative D4
	v4 = D4*v;
	% spectral derivative D4
	v = reshape(v,n,n);
	[D4x D4y] = df_2d([n n],4);
	v4x = ifft( D4x*fft(v) )/L^4;
	v4y = ifft( (D4y*fft(v')')' )'/L^4;
	v4 = reshape(v4_x+v4_y,n^2,1);
%}

	E	
	Err_e
	Err_v
	Err_est
end
figure(1); clf();
subplot(2,2,1);
loglog(N,abs(Err_e),'.-');
legend('location','southwest',leg_e)
xlabel('total number of grid points n');
ylabel('$\frac{ | \lambda_h {-} \lambda_* | }{|\lambda_*|}  $', 'interpreter', 'latex')
ylim([1e-2 1e0])
grid on
set(gca,'minorgrid','none');

subplot(2,2,2);
loglog(N(:,1),abs(Err_est),'.-');
legend('loation','southeast',leg_est{:})
xlabel('total number of grid points n');
ylim([1e-2 1e0])
grid on
set(gca,'minorgrid','none');

if (nargin() > 0 && printflag)
	preparePrint();
	print -depsc ../img/uniform-convergence.eps
	system('epstopdf ../img/uniform-convergence.eps');
end

subplot(2,2,3);
loglog(N,abs(Err_v),'.-');
legend(leg_v)
xlabel('total number of grid points n');
grid on
set(gca,'minorgrid','none');

figure(2);
subplot(2,2,1);
imagesc(reshape(-v_true,n,n))
colorbar()
subplot(2,2,2);
imagesc(reshape(v_fdm,n,n))
colorbar()
subplot(2,2,3);
imagesc(reshape(v_fem_1,n,n))
colorbar()


%loglog(M(:,2),Err_,'k'); %hold off
%loglog(M(:,1),abs(Err_est));
%legend({leg{:}, leg_est{:}} );
%grid on
%subplot(1,2,2);
%loglog(M(:,1),abs(Err_est));
%legend(leg_est);


end % function

