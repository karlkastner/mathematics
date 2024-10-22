% Wed 31 Jan 13:54:26 CET 2024
% Karl Kastner, Berlin
function test

	javaaddpath('./');
	rng(0);

	tab = table();

	n = 2^8*[1,1];
	L = n;
	dt = 1;
	nvar = 3;

	% diffusion coefficient
	% add 1 for dominance of diffusion
	d = dt*100*(1+rand(2,nvar));

	% advection coefficient
	ad = zeros(2,3,nvar);

	% reaction coefficient
	% add I for diagonal dominance
	aa = eye(nvar)+rand(nvar);
	a = {};
	% spatially varying reactiong coefficients
	sin_ = rand(nvar);
	afun = @(id,jd,n) ( aa(id,jd)*ones(n) ...
			    + sin(2*pi*id*(0:(n(1)-1))'/n(1)) ...
		  	    * sin(2*pi*jd*(0:(n(2)-1))/n(2)) ...
			  );
	for idx=1:nvar
		for jdx=1:nvar
			a{idx,jdx} = afun(idx,jdx,n);
		end
	end
	%d = dt*100*[1,1,100]; %(1:nvar);

	x0 = zeros(n(1),n(1),nvar);
	xi = zeros(n(1),n(1),nvar);
	b = zeros(n(1),n(2),nvar)+0.1*ones(n(1),n(1),nvar);
	for idx=1:nvar
		b(round(n(1)*idx/(nvar+1)),round(n(2)*idx/(nvar+1)),idx) = idx;
	end

	A = rad2d_assemble_matrix(a,ad,d,n,L,nvar);

	tic()
	%x_ref = A \ b(:);
	x_ref = b(:);
	tab.rt(1) = toc()
	x_ref = reshape(x_ref,[n,nvar]);

	%a0 = 1;

	mg_m = Multigrid();
	mg_m.init(a,ad,d,L,n,nvar);
	
	tic()
	[x_mg,resn]=mg_m.mldivide(b,x0);
	tab.rt(2) = toc()
	length(resn)
	resn(end)
	e.mg_m = rms(x_mg-x_ref,'all')
	max(abs(x_mg-x_ref),[],'all')

	%
	mg_2 = Multigrid();
	mg_2.init_fun(@(L,n,nvar,xi) rad2d_coeffs2diags2(d,ad,afun,L,n,nvar,xi),L,n,nvar,xi);
	tic()
	[x_2,resn]=mg_m.mldivide(b,x0);
	tab.rt(3) = toc()
	length(resn)
	resn(end)
	e.mg_2(1) = rms(x_2-x_ref,'all')
	e.mg_2(2) = rms(x_2-x_mg,'all')
	max(abs(x_2-x_ref),[],'all')

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


	if (0)
	figure(1)
	clf
	

	for idx=1:nvar
	subplot(nvar,4,1+(idx-1)*4)
	imagesc(x_ref(:,:,idx))
	axis equal;	
	axis tight	
	subplot(nvar,4,2+(idx-1)*4)
	imagesc(x_mg(:,:,idx))	
	axis equal;	
	axis tight	
	subplot(nvar,4,3+(idx-1)*4)
	imagesc(x_j(:,:,idx))	
	axis equal;	
	axis tight	
	subplot(nvar,4,4+(idx-1)*4)
	imagesc(x_gmres(:,:,idx))
	axis equal;
	axis tight	
	end
	end	
%	subplot(2,3,5)
%	imagesc(x_mg(:,:,2))	
%	subplot(2,3,3)
%	imagesc(x_gmres(:,:,1))	
%	subplot(2,3,6)
%	imagesc(x_gmres(:,:,2))	
	%imagesc(x_mg-x_ref);
	%colorbar

end % function test

