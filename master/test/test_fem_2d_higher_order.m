% Tue May  8 15:16:24 MSK 2012
% Karl Kästner, Berlin

function test_fem_2d_higher_order(printflag, N_)
	
	% number of grid points
	if (nargin() < 2 || isempty(N_))
		N_ = [2 4 8 16 32] + 1;
	end
	e_true = -2*pi^2;
	n_max = 1e4;
	lw = 2;
	
	% domain size
	L0 = [1 1];

	% measured values	
	E = [0];
	N = [0];
	Ta = [0];
	Te = [0];
	
	% constant grid loop
	for idx=1:length(N_)
		n = N_(idx)*[1 1];
	
		%	order, int-rule, legeng-entry
		C = {	1,     @int_2d_nc_3, 'p=1 Trapezoidal';
			1,  @int_2d_gauss_1, 'p=1 Gauss',
			2,  @int_2d_gauss_6, 'p=2 Gauss',
			3,  @int_2d_gauss_6, 'p=3 Gauss',
			4, @int_2d_gauss_12, 'p=4 Gauss',
			5, @int_2d_gauss_25, 'p=5 Gauss',
			6, @int_2d_gauss_33, 'p=5 Gauss',
		};
		int_exact = @int_2d_gauss_33;
	                   
		% get the grid
		[P T Bc] = mesh_2d_uniform(n,L0);
		mesh = Mesh_2d(P,T,Bc);
	
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
	
	log(Err_true(1:end-1,:)./Err_true(2:end,:))./log(N(1:end-1,:)./N(2:end,:))
	log(Err_est(1:end-1,:)./Err_est(2:end,:))./log(N(1:end-1,:)./N(2:end,:))
	
	%E = abs( (E - E_true)./E_true );
	Err_true
	%E - ones(size(E,1),1)*E(end,:)
	Err_est
	Err_est./Err_true
	Err_quad = E_ - e_true
	
	Tr = Ta+Te;
	
	figure(1); clf
	%subplot(2,2,1)
	subplot(1.5,1.5,1)
	loglog(N, Err_true,'.-','linewidth',lw); hold on
	loglog(N(1,:),Err_true(1,:),'k.--','linewidth',lw);
	loglog(N, Err_est,'.--','linewidth',lw); hold off
	xlabel('n : total number of grid points');
	ylabel('| \lambda_* - \lambda_h| / |\lambda_*|');
	ylim([1e-10 1e2]) 
	legend('Location', 'SouthWest', {C{:,3}, 'p increasing'}); 
	grid on;
	set(gca,'minorgrid','none');
	set(gca,'ytick',10.^(-10:2));

	if (nargin() > 0 && printflag)
		preparePrint();
		print -depsc ../img/fem-laplacian-2d-convergence.eps
		system('epstopdf ../img/fem-laplacian-2d-convergence.eps');
		return
	end
	
	subplot(2,2,3)
	loglog(N, abs(Err_quad),'-'); hold off
	grid on;
	set(gca,'minorgrid','none');
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
	ylim([1e-10 1e2]) 
	
	% adaptive grid refinement
	
end % function test_2d_higher_order

% TODO, NM into Mesh, promote into mesh
function [e err_est ta te n] = solve(P,T,Bc,f_promote,f_int,n_max)
	bcflag = 1;

	% construct 2nd order matrices
	tic()
	mesh = Mesh_2d(P,T,Bc);
	mesh.promote(f_promote);
	P = mesh.P;

	if (size(P,1) > n_max)
		e = NaN;
		err_est =NaN;
		ta =NaN;
		te = NaN;
		n = NaN;
		return;
	end

	mesh.prefetch();
	mesh.element_neighbours();

	% assemble discretisation matrices
	A = assemble_2d_dphi_dphi_java(mesh, [], f_int);
	B = assemble_2d_phi_phi_java(mesh, [], f_int);

	% apply boundary conditions
	[A B p__] = boundary_2d(A, B, mesh.Bc, bcflag);
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
		v_ = zeros(size(mesh.P,1),1); v_(p__,1) = v; v=v_;
	end
	
	% estimate the partial derivatives of degree p+1
	v_true = 2*sin(pi*P(:,1)).*sin(pi*P(:,2)); % why 2? v_ = v_/norm(v_);
if (f_promote < 0)
	dV  = mesh.dV(v, f_promote);
	% estimate the error
	obj = mesh.estimate_error(dV, 1, f_promote+1, 0);
	v_err = obj(1); err_est = obj(2); thresh = obj(3); nH = obj(4);

%	subplot(2,1,1)
%	display_2d(mesh, 0, [v(T(:,1),1).^2 v(T(:,2),1).^2 v(T(:,3),1).^2], [], 'EdgeColor', 'none');
%	caxis([0 20])
%	subplot(2,1,2)
%	display_2d(mesh, 0, [v_(T(:,1),1).^2 v_(T(:,2),1).^2 v_(T(:,3),1).^2], [], 'EdgeColor', 'none');
%	caxis([0 20])

%	[norm(v,'inf') max(nH) norm(double(dV),'inf')] %max(double(mesh.degen)) min(double(mesh.h_max))]
%	[max(nH) max(abs(dV))] %max(double(mesh.degen)) min(double(mesh.h_max))]
	% mark elements for refinement
	M = mark(v_err, thresh);
else
	err_est = NaN;
end

	n = size(A,1);
end

