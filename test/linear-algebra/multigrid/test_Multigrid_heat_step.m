% Mon  8 Jan 16:01:15 CET 2024
% v-cycle, stagnates
% full mg: hardly reduces number of iter
%
% -> (I - o*D^-1 A)
% convergence if 
%	0 < o < 2/max E(D^-1 A)
%	(1 + dt*2/dx^2)^-1*(a dt /dx^2, (1 - 2 a dt/dx^2),  a dt/dx^2)
%	E = 1 +/- 2 a dt/dx^2/(1 + a dt*2/dx^2)
%	  = 1 +/- 

	nd = 2;
	dt = 5e4;
%	dt = 1;
	reltol  = 1e-7;
	reltol  = 1e-5;
	maxiter = 1e4;
	o = 1;

	dx0 = 1
%	nn = [8,32,64,128,256,512,1024];
	nn = 128;
%nn = 32;
%	nn = [4,8,16,32,64,128];
%	nn = 2.^(2:12);
%	nn = 2.^12;
	t  = [];
	niter = [];
	L = max(nn)*dx0-1;
	L  = 1024;

	niter= [];
	rmse= [];
	rt = [];
	cond = [];
	oopt = [];
	for idx=1:length(nn)

	 n  = nn(idx);
	disp(n);
	 dx = L./n;
	 a0 = 1;
	 a  = 1;

	if (1 == nd)
		D2 = derivative_matrix_2_1d(n,L,2,'hdirichlet');
	else
		n = n*[1,1];
		%D2x, = derivative_matrix_2_1d(n,L,2,'hdirichlet');
		[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,L*[1,1],2,{'circular','circular'});
		D2=D2x+D2y;
	end
	I  = speye(prod(n));

%if (0)
%	 xi=(0:n-1)'/(n);
%	 b =sin(2*pi*xi);
%else
	if (1 == nd)
	 	b = zeros(n,1);
		b(n/2)=1;
	x0 = zeros(n,1);
	else
	 b = zeros(n);
	x0 = zeros(n);
		b(n(1)/2,n(2)/2)=1;
	end
%end



	A  = a0*I - (a*dt)*D2;
	D  = diag(A);
	iD = diag(sparse(1./D));
	R  = A-diag(sparse(D));
	sR = full(sum(abs(R),2));
	ll = iD(1)*(D(1)+sR(1)*[-1,+1]);
	oopt(idx,1) = 2./(ll(1)+ll(2));
	oopt(idx,2) = 2./ll(2);

if (0) % 1d only
	[min(D-sum(abs(R),2)),max(D+sum(abs(R),2))]
%	eest = eigs(iD*A,1,'lm')
	emax = 1 + 2*a*dt/dx^2/(1 + a*dt*2/dx^2)
end

	if (1)
	timer = tic();
	x_ref = A \ b(:);
	if (nd>1)
	x_ref = reshape(x_ref,n);
	end
%	[rcond(A), rcond(
%	cond(idx,1) = condest(iD*A);

	niter(idx,1) = 0;
	rt(idx,1) = toc(timer);
	rmse(idx,1) = 0;
	end
%	plot([b,x_ref])

if (0)
	for gamma=1:2
	timer = tic();
	if (1 == nd) 
		[x__,resn__]  = mg_heat_step_1d(a0,a,b,x0,dt,dx);
	else
		[x__,resn__]  = mg_heat_step_2d(a0,a,b,x0,dt,dx*[1,1],maxiter,reltol,gamma);
	end
	rt(idx,1+gamma)     = toc(timer);
	niter(idx,1+gamma)  = length(resn__);
	rmse(idx,1+gamma) = rms(flat(x__-x_ref));
	end
end

	% x_ = mg_heat_step(a0,a,b,x,n,dt,dx,x_ref,o);
%'molch'
	%if (1==nd)
	timer = tic();
	if (1==nd)
	[x_, resn] = mg_heat_1d_simple(a0,a,b,x0,n,dt,dx,reltol,maxiter);
	else
	mg = Multigrid();
	mg.init(a0,dt*a,L,n,1);
	
	[x_,resn]=mg.mldivide(b,x0);
	resn
rms(x_-x_ref,'all')
pause
	%[x_, resn] = mg_heat_2d_simple(a0,a,b,x0,n,dt,dx,reltol,maxiter);
	end
	rt(idx,3) = toc(timer);
	niter(idx,3) = length(resn);
	rmse(idx,3) = rms(flat(x_-x_ref));
	%end

%	niter(idx) = length(resn);
%	t(idx,2) = toc;
end
%	figure(1)
%	clf()
%	plot(resn)
%	hold on	
%	length(resn__)
%	loglog(nn,t)
	figure(1);
	clf
	subplot(2,2,1)
	plot(niter)
	subplot(2,2,2)
	loglog(nn,rt)
	subplot(2,2,3)
	plot(rmse)

	subplot(2,2,4)
	plot(oopt);
	if (nd == 1)
	%plot([x_ref,x_,x__])
	end

	figure(2)
	subplot(2,2,1)
	imagesc(x_ref)

	subplot(2,2,2)
	imagesc(x__)

	subplot(2,2,3)
%	imagesc(x_)
rmse(end,:)
if (0)
clf
	plot([x_ref,x_])
end

if (0)
pre  = 1; % Number of presmoothing iterations
post = 1; % Number of postsmoothing iterations
cycle = 1; % Type of multigrid cycle (1=V-cycle, 2=W-cycle, 3=F-cycle)
smooth = 1; % Smoother type (1=Jacobi, 2=Gauss-Seidel)
grids = 5; %log2(n)-2; % Max grids in cycle
maxit = 1e4; % Max iterations of solver
tol = 1e-07; % Tolerance of solver


[x_V,res_V,it_V] = multigrid(A,b,pre,post,cycle,smooth,grids,maxit,tol);
rms(x_V-x_ref)
rms(res_V)
rms(b)
it_V
	 plot([x_ref,x_,x_V])
end

% D^2 u = b

