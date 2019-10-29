% Sun Mar 25 14:07:34 MSK 2012
% Karl Kästner, Berlin

function test_radial_adaptive()
	path(path,'fem');
        path(path,'potential');
	m = 1;	% number of eigs to comp
	n = 4; % 3
	L0 = 100; %400;
	% initial grid including boundary points
	X  = L0*(0:n-1)'/(n-1);

	% initial shift
	s = -0.5;
	
	E = [];
	Err = [];

%	for idx=1:30
	idx=0;
	flag = 0;
%	while (idx<10)
	while(length(X) < 2.5e4 && (idx < 2 || abs(Err(idx,3)) < abs(Err(idx-1,3))) && 0 == flag)
		length(X)
		idx=idx+1;
		% analytic solution
		v_true = (sqrt(2)/sqrt(1)*exp(-X).*X); v_true = v_true/norm(v_true);
% fem
		[As B] = setup_fem(X);
		% find dominant eigenvalue
		[v_ e] = eigs(As-2*s*B,B,1,'SM'); v_=v_/norm(v_);
		E(idx,2) = 0.5*(e+2*s);
		% estimate the error
		Err(idx,4) = norm(v_(2:end-1).^2-v_true(2:end-1).^2)/norm(v_true.^2);
		Err(idx,5) = (-0.5 - E(idx,2))/0.5;
% fdm
		[As B D L_] = setup_fdm(X);
		% find dominant eigenvalue
		[v e flag] = eigs(As-s*B,B,1,'SM');
		E(idx,1) = e+s
		% undo similarity transform
		v = D*v;
		% adaptively refine the mesh
		X_=X;
		[X err v4 ] = fdm_adaptive_refinement({X}, length(X), v, L_, 2);
		X = X{1};
		% estimate the error
		v=v/norm(v);
		Err(idx,1) = norm(err);
		Err(idx,2) = norm(v.^2-v_true(2:end-1).^2)/norm(v_true.^2);
		Err(idx,3) = (-0.5 - E(idx,1))/0.5;
% plots
		figure(1)
		subplot(4,1,1)
		plot(X_,[[v_true.^2] [0; v.^2; 0] v_.^2])
		subplot(4,1,2)
		plot(X_(2:end-1),L_*v)
		subplot(4,1,3)
		plot(X_,v4)
		plot(X_(2:end-1),L_*L_*v)
%		semilogy(abs(fft(v4)))
		subplot(4,1,4)
		plot(X_(1:end-1),err)
%		subplot(5,1,5)
%		plot(X,1:length(X));
		subplot(4,1,1); hold on
		plot(X_,v_true.^2,'g');
		hold off
		N(idx) = length(X)
		norm(v)
		norm(v_)
	end % for idx
	figure(2);
	%loglog(N, abs(Err(:,1:3)),'.-')
	subplot(1.5,1.5,1)
	loglog(N, abs(Err(:,3)),'-ok','Linewidth',0.66,'Markersize',2.5,'Markerfacecolor','k')
	grid on; set(gca,'minorgrid','none')
	%xlabel('number of grid points n')
	xlabel('number of grid points n','interpreter','latex')
	ylabel(['$\left ( \sum_{i=1}^{' num2str(m) '} \left (\lambda_i {-} \lambda_i^{(n)} \right )^2 \right )^{1/2}$'],'interpreter','latex')
	legend('location','southwest','adaptive grid')
%	print -deps ../img/radial_adaptive.eps
%	subplot(3,3,7)
	figure(3)
	loglog(N,abs(Err(:,1:5)))
	legend('est','fdm v','fdm e','fem v', 'fem e')
%	semilogy(X(1:end-1),abs(diff(X)))
end % function test_radial_adaptive

function [As B D L_] = setup_fdm(X)
	[Ls D L_] = d_vargrid(X,2,2);
%	[L D] = laplacian_non_uniform(X); L_=D*L; D=sqrt(D); Ls = D*L*D;
	% potential
	V = diag(sparse(1./X(2:end-1)));
	% stiffness matrix
	As = -0.5*Ls - V;
	% enforce symmetry (lost due to rounding errors)
	As = 0.5*(As+As');
	% mass matrix
	B = speye(length(X)-2);
end

function [A B] = setup_fem(X)
	[P T] = mesh_1d_uniform(X);
	[A B] = fem_1d(P,T,@int_1d_cp,@potential);
end

function p = potential(x0)
	p = -2./sqrt(x0.^2);
end

