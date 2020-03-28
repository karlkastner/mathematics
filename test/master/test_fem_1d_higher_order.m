% Tue May  8 15:16:24 MSK 2012
% Karl Kästner, Berlin

function test_fem_2d_higher_order(N_)

% number of grid points
if (nargin() < 1 || isempty(N_))
	N_ = 2.^(1:3) +1; % 12
end
e_true = -pi^2;
n_max = 1e4;

% domain size
L0 = [1];

k_max = inf;

% constant grid loop
for idx=1:length(N_)
	n = N_(idx);
	k = 0;

	% get the grid
	[P T Bc] = mesh_1d_uniform(n,L0);
	mesh = Mesh_3d(P,T,Bc);

	k = k+1;
	f_int     = @int_1d_nc_2;
	f_promote = 1;
	[E(idx,k) Err_est(idx,k) Ta(idx,k) Te(idx,k) N(idx,k)] = solve(P, T, Bc, f_promote, f_int, n_max);
	leg{k} = 'linear trapezoidal';

	k = k+1;
	f_int = @int_1d_gauss_1;
	f_promote = 1;
	[E(idx,k) Err_est(idx,k) Ta(idx,k) Te(idx,k) N(idx,k)] = solve(P, T, Bc, f_promote, f_int, n_max);
	leg{k} = 'linear gauss 1';

	k = k+1;
	if (k > k_max) continue; end
	f_int = @int_1d_gauss_2;
	f_promote = 1;
	[E(idx,k) Err_est(idx,k) Ta(idx,k) Te(idx,k) N(idx,k)] = solve(P, T, Bc, f_promote, f_int, n_max);
	leg{k} = 'linear gauss 2';

	k = k+1;
	if (k > k_max) continue; end
	f_int = @int_1d_gauss_3;
	f_promote = 2;
	[E(idx,k) Err_est(idx,k) Ta(idx,k) Te(idx,k) N(idx,k)] = solve(P, T, Bc, f_promote, f_int, n_max);

	leg{k} = 'quadratic gauss';
	k = k+1;
	if (k > k_max) continue; end
	f_int = @int_1d_gauss_4;
	f_promote = 3;
	[E(idx,k) Err_est(idx,k) Ta(idx,k) Te(idx,k) N(idx,k)] = solve(P, T, Bc, f_promote, f_int, n_max);
	leg{k} = 'cubic gauss';
 
	k = k+1;
	if (k > k_max) continue; end
	f_int = @int_1d_gauss_5;
	f_promote = 4;
	[E(idx,k) Err_est(idx,k) Ta(idx,k) Te(idx,k) N(idx,k)] = solve(P, T, Bc, f_promote, f_int, n_max);
	leg{k} = 'quartic gauss';

%	e_true
	E
end % for idx

% error
Err_true = abs(E - e_true);

log(Err_true(1:end-1,:)./Err_true(2:end,:))./log(N(1:end-1,:)./N(2:end,:))
log(Err_est(1:end-1,:)./Err_est(2:end,:))./log(N(1:end-1,:)./N(2:end,:))

%E = abs( (E - E_true)./E_true );
%E - ones(size(E,1),1)*E(end,:)
%Err_est./Err_true
Err_true
Err_est
Tr = Ta+Te
leg

figure(1); clf
subplot(2,2,1)
loglog(N, Err_true); hold on
loglog(N, Err_est,'--'); hold off
xlabel('total number of grid points');
ylabel('relative error');
%legend('3pt tmp','3pt smp','3pt gauss (2)','3pt trapez','6pt gauss','10pt gauss', '15pt gauss', '21pt','3pt 1/3*(2*trapez + 1*smp)');
legend(leg); 
grid on; set(gca,'minorgrid','none');

%subplot(2,2,3)
%loglog(

%figure(2); clf
subplot(2,2,2)
loglog(N, Tr);
xlabel('number of grid points per axis');
ylabel('runtime [s]');
legend(leg); 
%legend('3pt tmp','3pt smp','3pt gauss','3pt trapez','6pt gauss','10pt gauss', '15pt gauss','21pt');
grid on; set(gca,'minorgrid','none');

% efficency
subplot(2,2,4); cla();
colormap('default');
c = colormap();
%figure(3); clf
subplot(2,2,4); cla();
l = length(leg);
for idx=1:l
	p = (idx-1)/l
	c = [p (1-p)  1];
	loglog(Tr(:,idx),Err_true(:,idx),'color',c); %c(idx,:));
	hold on;
end
hold off;
legend(leg); 
xlabel('runtime [s]');
ylabel('relative error');

% adaptive grid refinement

end % function test_2d_higher_order

% TODO, NM into Mesh, promote into mesh
function [e err_est ta te n] = solve(P,T,Bc,f_promote,f_int,n_max)
	bcflag = 1;

	% construct 2nd order matrices
	tic()
%mesh = Mesh_3d(P,T,Bc);
%mesh.element_neighbours();
%	Nm = mesh.Nm;
	switch (f_promote)
		case {1}
			% nothing to do
		case {2}
			[P T Bc] = promote_1d_2_3(P, T, Bc);
		case {3}
			[P T Bc] = promote_1d_2_4(P, T, Bc);
		case {4}
			[P T Bc] = promote_1d_2_5(P, T, Bc);
		case {5}
			[P T Bc] = promote_1d_2_6(P, T, Bc);
		otherwise
			st = dbstack();
			error(st.name,'Order of accuracy has to be an integer between 2 and 6');
	end
	if (size(P,1) > n_max)
		e = NaN;
		err_est =NaN;
		ta =NaN;
		te = NaN;
		n = NaN;
		return;
	end

	% assemble stiffness matrix
	A = assemble_1d_dphi_dphi(P, T, [], f_int);
	% assemble the mass matrix
	B = assemble_1d_phi_phi(P, T, [], f_int);

	% apply boundary conditions
	[A B p__] = boundary_1d(A, B, Bc, bcflag);
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
		v_ = zeros(size(P,1),1); v_(p__,1) = v; v=v_;
	end

	v_true = sin(pi*P(:,1));

	N = neighbour_1d(P,T)
	T
	[h_side C] = regularity_1d(P,T,Bc);
	[M err_est v_err] = mark_1d(P, T, v, N, h_side, C)

	n = size(A,1);
end

