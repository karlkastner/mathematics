% Tue May  8 15:16:24 MSK 2012
% Karl Kästner, Berlin

% TODO test error estimation and adaptive refinement
function test_fem_3d_higher_order(printflag, N_)

	% number of grid points
	if (nargin() < 2 || isempty(N_))
		N_ = 2.^(1:5) +1; % 12
	end
	e_true = -3*pi^2;
	n_max = 1e4;
	lw = 2;
	lim = [1e-6 1e2];
	
	% domain size
	L0 = [1];

	% measured values	
	E = [0];
	N = [0];
	Ta = [0];
	Te = [0];
	
	% constant grid loop
	for idx=1:length(N_)
		n = N_(idx);
	
		%	order, int-rule, legeng-entry
		C = {	1,     @int_3d_nc_4, 'p=1 Trapezoidal';
			1,  @int_3d_gauss_1, 'p=1 Gauss',
			2,  @int_3d_gauss_4, 'p=2 Gauss',
			3, @int_3d_gauss_11, 'p=3 Gauss',
			4, @int_3d_gauss_24, 'p=4 Gauss',
			5, @int_3d_gauss_45, 'p=5 Gauss',
		... %	6, @int_3d_gauss_45, 'p=6 Gauss'
		};
		int_exact = @int_3d_gauss_45;
	                   
		% get the grid
		[P T Bc] = mesh_3d_uniform([n n n],[L0 L0 L0]);
		mesh = Mesh_3d(P,T,Bc);
	
		for k=1:size(C,1)
			[E(idx,k) Err_est(idx,k) Ta(idx,k) Te(idx,k) N(idx,k)] = solve(P, T, Bc, C{k,1}, C{k,2}, n_max);
			if (~isequal(int_exact,C{k,2}))
				E_(idx,k) = solve(P, T, Bc, C{k,1}, int_exact, n_max);
			else
				E_(idx,k) = E(idx,k);
			end
		end
		E
	end % for idx
	
	% error
	Err_true = abs(E - e_true);
	
	%log(Err_true(1:end-1,:)./Err_true(2:end,:))./log(N(1:end-1,:)./N(2:end,:))
	%log(Err_est(1:end-1,:)./Err_est(2:end,:))./log(N(1:end-1,:)./N(2:end,:))
	
	%E = abs( (E - E_true)./E_true );
	%E - ones(size(E,1),1)*E(end,:)
	%Err_est./Err_true
	Err_true
	Err_quad = E_ - e_true
	%Err_est
	Tr = Ta+Te
	
	figure(1); clf
	%subplot(2,2,1)
	subplot(1.5,1.5,1)
	loglog(N, Err_true,'.-','linewidth',lw); hold on
	loglog(N(1,:),Err_true(1,:),'k.--','linewidth',lw);
	loglog(N, Err_est,'.--','linewidth',lw); hold off
	xlabel('n : total number of grid points');
	ylabel('| \lambda_* - \lambda_h| / |\lambda_*|');
	ylim(lim) 
	legend('Location', 'SouthWest', {C{:,3}, 'p increasing'}); 
	grid on;
	set(gca,'minorgrid','none');
	set(gca,'ytick',10.^(-10:2));

	if (nargin() > 0 && printflag)
		preparePrint();
		print -depsc ../img/fem-laplacian-3d-convergence.eps
		system('epstopdf ../img/fem-laplacian-3d-convergence.eps');
		return
	end
	
	subplot(2,2,3)
	loglog(N, abs(Err_quad),'-'); hold off
	grid on;
	set(gca,'minorgrid','none');
	
	%subplot(2,2,3)
	%loglog(
	
	%figure(2); clf
	subplot(2,2,2)
	loglog(N, Tr);
	xlabel('n : number of grid points per axis');
	ylabel('runtime [s]');
	legend('Location', 'NorthWest', C{:,3}); 
	grid on;
	set(gca,'minorgrid', 'none');
	
	% efficency
	subplot(2,2,4); cla();
	colormap('default');
	c = colormap();
	%figure(3); clf
	subplot(2,2,4); cla();
	l = size(C,1);
	for idx=1:l
		p = (idx-1)/l;
		c = [p (1-p)  1];
		loglog(Tr(:,idx), Err_true(:,idx), 'color', c);
		hold on;
	end
	hold off;
	legend(C{:,3});
	xlabel('runtime [s]');
	ylabel('| \lambda_* - \lambda_h| / |\lambda_*|');
	ylim(lim) 
	
end % function test_3d_higher_order

function [e err_est ta te n] = solve(P,T,Bc,f_promote,f_int,n_max)
	bcflag = 1;

	% construct 2nd order matrices
	tic()
	mesh = Mesh_3d(P,T,Bc);
	mesh.element_neighbours();
	mesh.promote(f_promote);

	if (mesh.np > n_max)
		e = NaN;
		err_est =NaN;
		ta =NaN;
		te = NaN;
		n = NaN;
		return;
	end
%[P T Bc] = get_mesh_arrays(mesh)
%size(T)
%size(diff(sort(T,2),2))
%T =sort(T,2);
%T(:,1:end-1) - T(:,2:end)
%sort(Bc,2)
%Bc(:,1:end-1) - Bc(:,2:end)
%f_promote

	mesh.prefetch();
%min(mesh.determinant)
	% assemble discretisation matrices
	A = assemble_3d_dphi_dphi_java(mesh, [], f_int);

	B = assemble_3d_phi_phi_java(mesh, [], f_int);

	% apply boundary conditions
	[A B p__] = boundary_3d(A, B, mesh.Bc, bcflag);
	ta = toc();

	% solve the eigenvalue problem
	tic();
	if (1 == size(A,1))
		[v e] = eig(full(A),full(B));
	else
		[v e] = eigs(A,B,1,'SM');
	end
	te = toc();
	if (bcflag)
		v_ = zeros(mesh.np,1); v_(p__,1) = v; v=v_;
	end

	[P T Bc] = get_mesh_arrays(mesh);
	v_true = sqrt(8)*sin(pi*P(:,1)).*sin(pi*P(:,2)).*sin(pi*P(:,3));
	err_est = NaN;
%{	
	% estimate the partial derivatives of degree p+1
	dV  = mesh.dV(v, f_promote);
	% estimate the error
	obj = mesh.estimate_error(dV, 1, f_promote+1, 0);
	v_err = obj(1); err_est = obj(2); thresh = obj(3); nH = obj(4);
	err_v_true = v_true - v;
%	[max(nH) norm(double(dV),'inf')/mesh.np; max(double(mesh.degen)) min(double(mesh.h_max))]
	% mark elements for refinement
%	M = mark(v_err, thresh);
%	err_est = 0;
	[ norm(v,'inf') max(nH) norm(double(dV),'inf') max(double(mesh.degen)) min(double(mesh.h_max))]
%}

	n = size(A,1);
end

