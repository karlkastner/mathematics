% Wed 31 Jan 13:54:26 CET 2024

	javaaddpath('./');

	tab = table();
	iter = [];
	rtt = [];

	n  = 2^10*[1,1];
	L  = n;
	dt = 1;
	nvar = 1;
	aa = [1,0.9,0.3;
              0.5,1,0.1;
              0.2,0.4,1];
	aa = aa(1:nvar,1:nvar);
	aa = eye(nvar);
	s  = 0.0;

	se = logspace(0,4,5+4);
	% note the number of iterations increase proportional to max(sex/sey,sey/sex)^0.5
	% which means that the runtime becomes infinite when there is only diffusion in one direction

	for edx=1:length(se)

	ex = 1; %/sqrt(se(edx));
	ey = 1/se(edx);
%	ex = se(edx); %/sqrt(se(edx));
%	ey = se(edx); %1/se(edx);
	
	aa = aa(1:nvar,1:nvar);
	a = {};
	for idx=1:nvar
		for jdx=1:nvar
			a{idx,jdx} = aa(idx,jdx)*ones(n) + s.*(rand(n)-0.5);
		end
	end
	%d = dt*100*(1:nvar);

	x0 = zeros(n(1),n(1),nvar);
	b  = zeros(n(1),n(2),nvar)+0.1*ones(n(1),n(1),nvar);
	for idx=1:nvar
		b(round(n(1)*idx/(nvar+1)),round(n(2)*idx/(nvar+1)),idx) = idx;
	end

	A = assemble_rad_nvar(a,ex,ey,n,L);

	tic()
	x_ref = A \ b(:);
	%x_ref = b(:);
	rtt(edx,1)=toc();
	x_ref = reshape(x_ref,[n,nvar]);

	mg_m = Multigrid();
	mg_m.init(a,[ex;ey],L,n,nvar);
	
	tic();
	mg_m.maxiter = 1e3;
	[x_mgm,tab.flag(2),tab.rens(2),iter(edx,2)] = mg_m.mldivide(b,x0);
	x_m = b;
	%resn = NaN
	rtt(edx,2) = toc();
	%length(resn)
	%resn(end)
	tab.e(2) = rms(x_mgm-x_ref,'all')
	max(abs(x_mgm-x_ref),[],'all')

	%iter(edx,1) = tab.iter(2);

	mgmfun = @(x) flat(mg_m.cycle1(reshape(x,n(1),n(2)),0.*b));

	%reshape(x,[],obj.nvar)';
	%x = obj.aux.mg_j.mldivide(x,x0);
	%x=flat(x');

	% obj.aux.mg_j.set_a(obj.aux.aa,q*dt*rvec(obj.p.ex));
	% TODO preocompute initial vector or avoid it altogether
	%prec = {@mgmfun};
	%otherwise
	%		prec = {};
	%end
	%dz = bicgstabl(A,g,obj.opt.inner2_tol,maxiter,prec{:});

	tic();
	[x_,flag_,relres_,iter_] = gmres(A,b(:),[],tol_bicg,maxiter_bicg);
	rtt(edx,3)  = toc();
	iter(edx,3) = iter_(2);

	tic();
	[x_,flag_,relres_,iter(edx,4)] = bicgstabl(A,b(:),tol_bicg,maxiter_bicg);
	rtt(edx,4) = toc();

	maxiter_bicg = 1000;
	tol_bicg     = sqrt(eps);
	%GMRES(A,B,RESTART,TOL,MAXIT,M)
	tic();
	[x_,flag_,relres_,iter_] = gmres(A,b(:),[],tol_bicg,maxiter_bicg,mgmfun);
	rtt(edx,5) = toc();
	iter(edx,5)=iter_(2);

	tic();
	[x_,flag_,relres_,iter(edx,6)] = bicgstabl(A,b(:),tol_bicg,maxiter_bicg,mgmfun);
	rtt(edx,6) = toc();

	iter
	rtt
%end
	end % for edx

	figure(2);
	clf();
	subplot(2,2,1);
	loglog(se,iter);
	xlabel('max(ex,ey)/min(ex,ey)');
	ylabel('iter/iter(1)');

	subplot(2,2,2);
	loglog(se,rtt);
	legend('direct','mg','gmres','bicgstab','mg-gmres','mg-bicgstab');

if (0)
	disp('java')
	mg_j = javaObject('Multigrid_java');
	tic();
	aa_ = zeros(nvar,nvar,n(1)*n(2));
	for idx=1:nvar
	for jdx=1:nvar
		aa_(idx,jdx,:) = flat(a{idx,jdx});
	end
	end
	mg_j.init(aa_,d,L,n);
	%mg_j.init(aa,d,L,n);
	[x_j]=mg_j.mldivide(reshape(b,[],nvar)',reshape(x0,[],nvar)');
	tab.rt(3) = toc()
	x_j = reshape(x_j',n(1),n(2),nvar);
	iter = mg_j.iter
	resn = mg_j.resn
	e.mg_j = rms(x_j-x_ref,'all')


	tic()
	[Lp,Up] = ilu(A);
	restart = 10;
	%x_gmres = gmres(A,b(:),restart,[],1000,Lp,Up);
	x_gmres = b(:);
	tab.rt(4) = toc()	
	x_gmres = reshape(x_gmres,[n,nvar]);
	e.gmres = rms(x_gmres-x_ref,'all')
	
end	

	figure(1)
	clf();
	

	for idx=1:nvar
	subplot(nvar,4,1+(idx-1)*4)
	imagesc(x_ref(:,:,idx))
	axis equal;	
	axis tight
	title('x_ref')
	
	subplot(nvar,4,2+(idx-1)*4)
	imagesc(x_mgm(:,:,idx))	
	axis equal;	
	axis tight	
	title('x_mgm')
	
	subplot(nvar,4,3+(idx-1)*4)
%	imagesc(x_j(:,:,idx))	
	axis equal;	
	axis tight	
	subplot(nvar,4,4+(idx-1)*4)
%	imagesc(x_gmres(:,:,idx))
	axis equal;
	axis tight	
	end	
%	subplot(2,3,5)
%	imagesc(x_mg(:,:,2))	
%	subplot(2,3,3)
%	imagesc(x_gmres(:,:,1))	
%	subplot(2,3,6)
%	imagesc(x_gmres(:,:,2))	
	%imagesc(x_mg-x_ref);
	%colorbar
	
