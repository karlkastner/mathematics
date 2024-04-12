% Wed 31 Jan 13:54:26 CET 2024

	javaaddpath('/home/pia/phd/src/lib/mathematics/linear-algebra/multigrid/');

	tab = table();
	iter = [];

	n  = 2^10*[1,1];
	L  = n;
	dt = 10;
	nvar = 1;
	aa = [1,0.9,0.3;
              0.5,1,0.1;
              0.2,0.4,1];
	aa = aa(1:nvar,1:nvar);
	aa = eye(nvar);
	s  = 0.0;
	ad = zeros(2,3);
	ad(1,2:3) = [-1,1]/2;
	ad(2,2:3) = 2*[-1,1];
%	ad(1,1:2) = [1,-1];
	exy = [1/2;2];
	%exy = [0;1];
	%exy = [2;1/2];
	%ad(1,2) = -0.1;
	%ad(1,1) = -0.2;

	se = logspace(0,2,5+4);
	% note the number of iterations increase proportional to max(sex/sey,sey/sex)^0.5
	% which means that the runtime becomes infinite when there is only diffusion in one direction

	for edx=1:1 %length(se)

	
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

	A = assemble_rad_nvar(a,dt*ad,dt*exy,n,L,nvar);

	tic()
	x_ref = A \ b(:);
	%x_ref = b(:);
	tab.rt(1) = toc()
	x_ref = reshape(x_ref,[n,nvar]);

	mg_m = Multigrid();
	mg_m.init(a,dt*ad,dt*exy,L,n,nvar);
	
	tic()
	mg_m.maxiter = 1e4;
	[x_mgm,tab.flag(2),tab.relres(2),tab.iter(2)] = mg_m.mldivide(b,x0);
	x_m = b;
	%resn = NaN
	tab.rt(2) = toc()
	%length(resn)
	%resn(end)
	tab.e(2) = rms(x_mgm-x_ref,'all')
	%max(abs(x_mgm-x_ref),[],'all')
	rms(x_ref-flipud(x_mgm),'all')
	resn0 = rms(b(:));
	resn = rms(A*flat(x_mgm)-b(:));
	resn/resn0
	tab
	
	Ax=reshape(A*flat(x_mgm),n);
	mg_m.s(1).b = 0;
	mg_m.s(1).x = x_mgm;
	mg_m.resfun(1);
	Ax_ = mg_m.s(1).res;
	rms(Ax-Ax_,'all')
%pause
	%iter(edx,1) = tab.iter(2);

	mgmfun = @(x) flat(mg_m.cycle1(b,reshape(x,n(1),n(2))')');

	x = b;
	for jdx=1:16
		x = mg_m.cycle1(b,x);	
	end
	rms(x_ref-x,'all')

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

if (0)
	maxiter_bicg = 100;
	tol_bicg     = sqrt(eps);
	[x_,flag_,relres_,iter(edx,2)] = bicgstabl(A,b(:),tol_bicg,maxiter_bicg,mgmfun);
	iter
end
if (1)
	%disp('java')
	mg_j = javaObject('Multigrid_java');
	tic();
	aa_ = zeros(nvar,nvar,n(1)*n(2));
	for idx=1:nvar
	    for jdx=1:nvar
		aa_(idx,jdx,:) = flat(a{idx,jdx});
	    end
	end
	aa_  = 1;
	mg_j.init(aa_,dt*ad,dt*exy,L,n);
	%mg_j.init(aa,d,L,n);
	tic();
	x_j = mg_j.mldivide(reshape(b,[],nvar)',reshape(x0,[],nvar)');
	%tab.rt(3) = toc()
	x_j = reshape(x_j',n(1),n(2),nvar);
	tab.rt(3) = toc();
	tab.iter(3) = mg_j.iter
%	resn = mg_j.resn
%	e.mg_j = rms(x_j-x_ref,'all')
	tab.e(3) = rms(x_j-x_ref,'all')
end

%end
	end % for edx

if(0)
	figure(2);
	clf();
	loglog(se,iter);
	xlabel('max(ex,ey)/min(ex,ey)');
	ylabel('iter/iter(1)');
end
if (0)

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
	

	nvar = 2;
	idx =1;
%	for idx=1:nvar
	subplot(nvar,4,1+(idx-1)*4)
	imagesc(x_ref(:,:,idx))
	axis equal;	
	axis tight
	title('x_ref')
	%imagesc(x_ref(:,:,idx)-flipud(x_mgm(:,:,idx)))
	colorbar	
	
	subplot(nvar,4,2+(idx-1)*4)
	imagesc(x_mgm(:,:,idx))	
	axis equal;	
	axis tight	
	title('x_mgm')
	%imagesc(x_ref(:,:,idx)-flipud(x_mgm(:,:,idx)))
	colorbar	
	
	subplot(nvar,4,3+(idx-1)*4)
	imagesc(x_j);
	%ref(:,:,idx)-(x_mgm(:,:,idx)))
	colorbar	
%	imagesc(x_j(:,:,idx))	
	axis equal;	
	axis tight	
	subplot(nvar,4,4+(idx-1)*4)
%	imagesc(x_gmres(:,:,idx))
	axis equal;
	axis tight

%	end	
	subplot(2,4,5)
	imagesc(Ax);

	subplot(2,4,6)
	imagesc(Ax_);


%	subplot(2,3,5)
%	imagesc(x_mg(:,:,2))	
%	subplot(2,3,3)
%	imagesc(x_gmres(:,:,1))	
%	subplot(2,3,6)
%	imagesc(x_gmres(:,:,2))	
	%imagesc(x_mg-x_ref);
	%colorbar
	
