% Sat Sep 10 14:45:25 MSD 2011
% Karl Kästner, Berlin
%
% compute energy levels of a confined hydrogen atom
% eigenvalues of the discretised schrödinger equation
%
% function energy_level(dimension, n_grid_0, n_refinement)
% dimension    : { 1, 2, 3 }
% n_grid_0     : number of grid points in one dimension to start with
% n_refinement : number of refinement steps
% k_refinement : factor to increase grid points per iteration (can be % irrational)
% lim          : maximum rank of hessenberg-matrix
function energy_level(dimension, n_grid_0, n_refinement, k_refinement, lim)

% number of grid points (must be even, as hydrogen atom is located in centre)
n_grid_0 = 2*floor(n_grid_0/2);
n_grid_d = n_grid_0;

% colours of plot
clf
hsv = [ (1:n_refinement+1)'/(n_refinement+1), ones(n_refinement+1,1), ones(n_refinement+1,1) ];
colour = hsv2rgb(hsv);
% legend entries
l = {};

for ndx=1:n_refinement+1
	n_grid = round(n_grid_d);
	n_grid = 2*floor(n_grid/2);

	[A h] = hydrogen_boxed(n_grid, dimensions);
	l{end+1} = sprintf('h = %e',h);

	% slow !
	% disp(sprintf('condest: L %e V %e A %e\n', condest(L), condest(V), condest(A)));
	
	% impose boundary conditions
	% nothing to do for dirichlet conditions

	% simplify problem by arnoldi method
	tic();
	if (size(A,1) > lim)
		A = arnoldi(A, lim+1, 1);
		A = A(1:end-1,1:end);
		% slow : test loss of orthogonality
		% rank(full(A))
	end
	t1 = toc();

	% compute eigenvalues = energy niveau
	tic()
	[V E F] = eigs(A, n_grid_0, 'SM');
	t2 = toc();
	disp([num2str(n_grid) ' ' num2str(h) ' ' num2str(size(A,1)) ' ' num2str(t1) ' ' num2str(t2)]);
	% remove spurious eigenvalues
	E = diag(E);
	E = unique(single(E));
	E = sort(E);
	E = E(1:min(n_grid_0,length(E)));

	if (0 ~= F)
		'no convergence of eigs'
		return
	end

	% plot result
	subplot(1,2,1)
	plot(1:min(n_grid_0,length(E)), E, '.','Markersize',24,'Color',colour(ndx,:));
	hold on

	% prepare refinement in next step
	n_grid_d = k_refinement*n_grid;
end % for ndx

title('Eigenvalues of \Delta\psi + (\alpha + 2/r)\psi = 0');
ylabel('value of eigenvalue')
xlabel('smallest eigenvalues, sorted ascending');
legend('Location','Best',l);
grid on
set(gca,'Minorgrid','none')

E

% plot eigenvectors
subplot(1,2,2)
plot(V)
grid on
%xlim([0 1])

end % function energy_level

