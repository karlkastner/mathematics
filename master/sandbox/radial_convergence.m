% 2012 Jan 22 00:40 (MSK)
% Karl Kästner, Berlin

% Radial Schrödinger:
% 1/2 d/dr ( r^2 d/dr R ) + r R = -E r^2 R
% Substitution: y = rR, y(0)=0, y(r->inf) = 0
% 1/2 y'' + 1/r y = -E y

function [EE F] = radial_convergence()
	figure(1); clf
	figure(2); clf
	path(path,'./derive')
	
	clf
	N  =  2.^(2:14);
	L0 =  400;
%	m=10;
	m=1;
	s = -0.5;
	f = { @setup_constgrid, @setup_vargrid }
%	f = { @setup_vargrid }
	
	P = [];
	
	for mdx=1:length(f)
		for ndx=1:length(N);
			n = N(ndx)
			[A B] = feval(f{mdx},n,L0);
			I = speye(size(A));
			E(1:min(n,m),ndx) = sort(eigs(A-s*I,[],min(n,m),'SM'))+s;
		end % for ndx
		
		%	todo undo R = 1/r y and also undo symmetry transform
		%	[v e] = eigs(A,B,min(n,n_c),'SM');
		%	[e idx] = sort(diag(e));
		%	v = v(:,idx);
		%	X = 400*(1:N(end))/N(end);
		%	plot(X,v(:,1).*diag(V))
		%	xlim([0 10])
		
		% richardson interpolation
		for ndx=1:length(N)
			c = derive_richardson(min(ndx,100))';
			F(:,ndx) = E(:,ndx+1-min(100,ndx):ndx)*c;
		end % for ndx
		
		G = -0.5./(1:m)'.^2*ones(1,length(N));
		
		E_err = E(1:m,:) - G;
		F_err = F(1:m,:) - G;
		
		P(:,1+2*(mdx-1)) = sqrt(sum(E_err.^2,1))';
		P(:,2+2*(mdx-1)) = sqrt(sum(F_err.^2,1))';
		figure(1);
		subplot(2,2,mdx)
		loglog(N,P(:,1+2*(mdx-1)),'k-','Linewidth',2); hold on
		loglog(N,P(:,2+2*(mdx-1)),'k--','Linewidth',2); hold on
		grid on
		legend('2nd order FDM', 'extrapolated')  
		if (1 == mdx)
			title('Radial Schrödinger Equation - Constant Grid')
		else
			title('Radial Schrödinger Equation - Variable Grid')
		end % if
		xlabel('number of grid points n')
		ylabel('residual norm of first 10 eigenvalues')
		set(gca,'minorgrid','none')
		set(gca,'xtick',2.^(0:16))
		set(gca,'ytick',10.^(-15:0))
		xlim([N(1) N(end)])
		ylim([1e-10 1e0])
		EE(1:size(E,1),1:size(E,2),mdx)=E
	end % for mdx
	
	figure(2)
	subplot(1.5,1.5,1)
	set(gca,'nextplot','replacechildren')
	set(gca,'LineStyleOrder',{'-o','--o','-^','--^'})
	set(gca,'ColorOrder',[0 0 0])

	loglog(N, P, 'Linewidth',0.66,'Markerfacecolor','k','Markersize',2.5) % 3
	set(gca,'yscale','log','xscale','log')
%	title('Radial Schr\"odinger Equation','interpreter','latex')
	legend('Location','SouthWest','uniform grid', 'u.g. extrapolated', 'varying grid', 'v.g. extrapolated','interpreter','latex')
	%xlabel('number of grid points n')
	xlabel('number of grid points n','interpreter','latex')
	%ylabel('$\left ( \sum_{i=1}^{10} \left (\lambda_i {-} \lambda_i^{(n)} \right )^2 \right )^{1/2}$','interpreter','latex')
	ylabel(['$\left ( \sum_{i=1}^{' num2str(m) '} \left (\lambda_i {-} \lambda_i^{(n)} \right )^2 \right )^{1/2}$'],'interpreter','latex')
	grid on
	set(gca,'minorgrid','none')
	%set(gca,'ytick',10.^(-10:1:0))
	set(gca,'xtick',10.^(0:1:10)) 
	set(gca,'box','on')
	%xlim([N(1) N(end)])
	ylim([1e-10 1e0])
	
	print -depsc ../img/radial-convergence-richardson.eps
end % radial_convergence

function [A B V] = setup_constgrid(n,L0)
	% laplacian
	L = (n+1)^2/L0^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	% potential
	V = diag(sparse((n+1)./(L0*(1:n))));
	% stiffness matrix
	A = -0.5*L - V;
	% mass matrix
	B = [];
end % setup_constgrid

function [As B] = setup_vargrid(n, L0)
	% variable grid
%	X = 0.5*(xgrid(n+2, 2, -L0, L0)' + L0);
	e = 5/L0;%0.5;
	h = L0/(n+1);
	X = h*(0:n+1); % - L0/2;
	X = (X(end)/(e*(exp(e*X(end))-1))) * e * sign(X).*(exp(e*abs(X)) - 1);
	X = X';
	% laplacian
%	Ls = laplacian_non_uniform(X);
	[L D1] = laplacian_non_uniform(X);
	Ls = sqrt(D1)*L*sqrt(D1);
	% potential
	V = diag(sparse(1./X(2:end-1)));
	% stiffness matrix
	As = -0.5*Ls - V;
	% enforce symmetrie
	As = 0.5*(As+As');
	% mass matrix
	B = [];
end % setup_vargrid

