% Wed 31 Jan 13:54:26 CET 2024

% benchmark: s for 1024^2
% n 	direct	mgm	mgj	gmres
% 1024	
%	N/A	16	14	10*17+4		gmres (10) with ilu
%	112.9	17.4	2.37	64.0		with side diagonals/cross reaction, 24g for direct
%
% 1024	diagonal only
% 	N/A 		21		18		10*14+2
% 	16.058854	23.540634	2.345259	53.680582
%	13.5		15.99		2.48
%	N/A		16		14
%	14.6		15.7		2.5		50.0		
% 512	
%	N/A	17	14	10*17 + 6
%	2.76	0.53	0.13	29.4
%
% 512	diag
%		14	17	10*15+8
%	3.4	2.8	0.6	15.1
%
% 512 random coeff, s = 0.1
%	N/A	19	17	10*17+9
%	16.4	3.5	0.88	18.2
% 	17.3	3.2	0.6	17.5
%	16.5	2.9	0.8	16.4	
%
% upscaling:
% t1		   0.6s
% single precision 1/2
% time steps	   5e5/1e3
% newton steps	   3
%		   1080s ~ 20 min vs 30 h vs 

	javaaddpath('./');

	tab = table();

	n = 2^10*[1,1];
	L = n;
	dt = 1;
	%a0 = 1 + dt*(rand(n)-0.5);
	%a0(:,2) = 1 + dt*(rand(n)-0.5);
	nvar = 3;
	aa = [1,0.9,0.3;
              0.5,1,0.1;
              0.2,0.4,1];
	aa = eye(3);
	%s = 0;
	s = 0.1;
	%s = rand(nvar,nvar);
	aa = aa(1:nvar,1:nvar);
	%aa = eye(nvar);
	a = {};
	for idx=1:nvar
	for jdx=1:nvar
		a{idx,jdx} = aa(idx,jdx)*ones(n) + s.*(rand(n)-0.5);
	end
	end
	%d = dt*100*[1,1,100]; %(1:nvar);
	d = dt*100*(1:nvar);

	x0 = zeros(n(1),n(1),nvar);
	%b  = zeros(n(1),n(1),2);
	%b(n/2,n/2) = 1;
	%b = ones(n);
	b = zeros(n(1),n(2),nvar)+0.1*ones(n(1),n(1),nvar);
	for idx=1:nvar
		b(round(n(1)*idx/(nvar+1)),round(n(2)*idx/(nvar+1)),idx) = idx;
	end
	%b(round(n(1)*2/3),round(n(2)*2/3),2) = 2;

	A = assemble_rad_nvar(a,d,n,L);

	tic()
	%x_ref = A \ b(:);
	x_ref = b(:);
	tab.rt(1) = toc()
	x_ref = reshape(x_ref,[n,nvar]);

	%a0 = 1;

	mg_m = Multigrid();
	mg_m.init(a,d,L,n,nvar);
	
	tic()
	%[x_mg,resn]=mg_m.mldivide(b,x0);
	%[x_mg,resn]=mg_m.mldivide(b,x0);
	x_mg = b; resn=NaN;
	tab.rt(2) = toc()
	length(resn)
	resn(end)
	e.mg_m = rms(x_mg-x_ref,'all')
	max(abs(x_mg-x_ref),[],'all')

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
%	subplot(2,3,5)
%	imagesc(x_mg(:,:,2))	
%	subplot(2,3,3)
%	imagesc(x_gmres(:,:,1))	
%	subplot(2,3,6)
%	imagesc(x_gmres(:,:,2))	
	%imagesc(x_mg-x_ref);
	%colorbar
	
