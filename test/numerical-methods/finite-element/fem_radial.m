% 2012 Apr 30 15:47 (MSK)
% Karl Kästner, Berlin


function fem_radial(printflag)
	
	N = 2.^(1:10)
	
	s=-1; %0.25;
	S=s;
	
	%% adaptive grid
	L = 400;
	P = [0; L/2; L];
	T = [1 2; 2 3];
	midx = 1;
	madx = 3;
	Bc = [1; 3];
	
	idx=1;
	while (length(P) < 1e4)
		% discretisation matrices
	%	[A B] = fem_radial_int(P, T, @int_1d_gauss_2);
	%	[A B] = fem_radial_int(P, T, @int_1d_cp);
	%	[A B] = fem_radial_int(P, T, @int_1d_corner);
		[sP sdx] = sort(P);
		[A B] = fem_radial_exact(sP);
	
		% apply bc
		% only for numerical integration
	%	a=A(midx,midx); a=1; A(midx,:) = 0; A(midx,1)=-a*1e12; B(midx,:) = 0; B(midx,midx) = a;
	%	a=A(madx,madx); a=1; A(madx,:) = 0; A(madx,1)=-a*1e12; B(madx,:) = 0; B(madx,madx) = a;
	
		% eigenvalues and vectors
		if (length(A) == 1)
		        [v e] = eig(full(A-s*B),full(B));
		else
		        [v e] = eigs(A-s*B,B,1,'SM');
		end
		% only for exact integration
		v = [0; v; 0];
		v(sdx) = v;
	
		E(idx,1) = e+s;
		M(idx,1) = size(T,1);
		v_true = P.*exp(-P);
		v_true = v_true/norm(v_true);
		v = v/norm(v);
		vErr(idx,1) = norm(v.^2 - v_true.^2);
	
		% adaptively refine the mesh
		Nm = neighbour_1d(P, T);
		[h_side C] = regularity_1d(P,T,Bc);
		[M_ nerr v_err] = mark_1d(P, T, v(:,1), Nm, h_side, C);
		[P T] = refine_1d(P, T, M_);
		%[P T nerr dv C ve] = fem_adapt_1d(P,T,v);
		Nerr(idx,1) = nerr;
		M(idx)
		idx=idx+1;
	end
	eErr = abs(E + 1);
	figure(1);
	set(gcf,'DefaultAxesColorOrder',[0 0 0]);
	set(gcf,'DefaultAxesLineStyleOrder',{'*','-^','-o','-s','+','-v','--','-'});
	subplot(1.5,1.5,1)
	loglog(M(2:end),eErr(2:end),'-ok'); hold on
	loglog(M(2:end),vErr(2:end),'-vk');
	loglog(M(2:end),Nerr(2:end),'-k'); hold off
	legend('location','southwest','eigenvalue','eigenvector','estimate');
	grid on
	set(gca,'xtick',10.^(0:5));
	set(gca,'ytick',10.^(-10:0));
	set(gca,'minorgrid','none')
	xlabel('number of grid points m')

	set(gca,'ColorOrder',[0 0 0]);
	set(gca,'LineStyleOrder',{'-*','-^','-o','-s','+','-v','--','-'});
	if (nargin() > 0 && printflag)
		preparePrint();
		print -deps ../img/fem-convergence-1d-adaptive.eps
		system('epstopdf ../img/fem-convergence-1d-adaptive.eps');
	end
	
	%% constant grid
	L = 10;
	E = [];
	for idx=1:length(N)
		n = N(idx);
		X = L*(1:n)'/(n+1);
		h = L/(n+1);
		v_true = X.*exp(-X);
		v_true = v_true / norm(v_true);
	
		F={@int_1d_gauss_1, @int_1d_nc_2, @int_1d_gauss_2};
		for fdx=1:length(F)
			[A B] = fem_radial_int([0; X; L], [(1:n+1)' (2:n+2)'], F{fdx});
			% BC
			%[A B] = boundary(A, B, Bc);
			A = A(2:end-1,2:end-1);
			B = B(2:end-1,2:end-1);
			% compute eigenvalues
			[v e] = eigs(A-s*B,B,1,'SM');
			v=v/norm(v);
			E(idx,1+fdx) = e+s;
			Verr(idx,1+fdx) = norm(v.^2-v_true.^2);
		end
	
		[A B] = fem_radial_exact([0; X; L]);
		% compute eigenvalues
		[v e] = eigs(A-s*B,B,1,'SM');
		v=v/norm(v);
		E(idx,1) = e+s;
		Verr(idx,1) = norm(v.^2-v_true.^2);
		E
	end
	Err = abs(E-S); %abs(E/pi^2-1)
	
	figure(2);
	set(gcf,'DefaultAxesColorOrder',[0 0 0]);
	set(gcf,'DefaultAxesLineStyleOrder',{'*','-^','-o','-s','+','-v','--','-'});
	subplot(2,2,1);
	loglog(N,Err);
	xlabel('number of grid points m');
	ylabel('|\lambda_* \lambda|/|\lambda_*|');
	title('First eigenvalue');
	grid on;
	set(gca,'minorgrid','none');
	legend('location','southwest','exact','midpoint','trapezoidal','2-point Gauss');
	xlim([1 N(end)]);
	ylim([1e-6 1]);
	set(gca,'xtick',10.^(0:3));
	
	subplot(2,2,2);
	loglog(N,Verr);
	xlabel('number of grid points m');
	ylabel('||v_* - v||_2');
	
	title('First eigenvector');
	grid on;
	set(gca,'minorgrid','none');
	legend('location','southwest','exact','midpoint','trapezoidal','2-point Gauss')
	%plot(X,[v_true.^2 - v.^2])
	xlim([1 N(end)]);
	ylim([1e-6 1]);
	set(gca,'xtick',10.^(0:3));
	
	if (nargin() > 0 && printflag)
		preparePrint();
		print -deps ../img/fem-radial-convergence.eps
		system('epstopdf ../img/fem-radial-convergence.eps')
	end
end % fem_radial

function [A B] = fem_radial_int(P,T,f_int)
	A = assemble_1d_dphi_dphi(P, T, [], f_int);
	V = assemble_1d_phi_phi(P, T, @potential_coulomb, f_int);
	B = assemble_1d_phi_phi(P, T, [], f_int);
	A = -A + 2*V;
end

function [A B] = fem_radial_exact(X)
	n = length(X)-2;
	xl = [X(1:end-2)];
	xc = X(2:end-1);
	xr = [X(3:end)];
%	h = X(2)-X(1); % todo, vargrid
%	h = X(2)-X(1); % todo, fix for variable grid-width
%	A = spdiags(ones(n,1)*1/h*[-1 2 -1],-1:1,n,n);
%	B = spdiags(ones(n,1)*h/6*[1 4 1],-1:1,n,n);
	% stiffness - Laplacian
	al = -1./(xc - xl);
	ar = -1./(xr - xc);
	ac = -(al+ar);
	A = spdiags([ [al(2:end); NaN] ac [NaN; ar(1:end-1)]], -1:1, n, n);
	% mass
	bl = 1/6*(xc - xl);
	br = 1/6*(xr - xc);
	bc = 1/3*(xr - xl);
	B = spdiags([ [bl(2:end); NaN] bc [NaN; br(1:end-1)]], -1:1, n, n);

	% stiffness - Potential
	Vl = -0.5*(xl+xc)./(xc-xl) + xc.*xl.*log(xc./xl)./(xc-xl).^2;
	Vr = -0.5*(xc+xr)./(xr-xc) + xr.*xc.*log(xr./xc)./(xr-xc).^2;
	Vc = -xl./(xc-xl) + (xl.^2.*log(xc./xl))./(xc-xl).^2 ...
	     -xr./(xr-xc) + (xr.^2.*log(xr./xc))./(xr-xc).^2;
	Vc(1,1) = -xr(1)/(xc(1)-xl(1)) + (xr(1).^2.*log(xr(1)./xc(1)))./(xr(1)-xc(1)).^2;
	V = spdiags([	[Vl(2:end); NaN], -Vc, [NaN; Vr(1:end-1)] ], -1:1, n, n );
%	full(A)
%	full(B)
%	full(V)
%	pause
	A = A + 2*V;
end

