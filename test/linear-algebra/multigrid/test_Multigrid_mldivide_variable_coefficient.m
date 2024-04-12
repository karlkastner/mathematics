% 2024-01-31 11:09:32.136836006 +0100
	nvar = 3;
	n = 2^7*[1,1];
	L = n;
	dt = 100;
	s = 0.1;
	d = 1;
	x = linspace(0,1,n(1));
	y = linspace(0,1,n(2))';
	%a = cos(pi*x+pi/2).*cos(pi*y+pi/2);
	d = (1:nvar);

	%a0 = 1 + 0.*dt*(rand(n)-0.5);
	%a0 = 1 + (a-0.5);
	%b = repmat(rvec(b(:)),1,nvar);
	%x0 = repmat(rvec(x0(:)),1,nvar);
	x0 = zeros(nvar*prod(n),1);
	b  = zeros(n(1),n(2),nvar);
	for idx=1:nvar
		b(round(idx*n/(nvar+1)),round(idx*n/(nvar+1))) = idx;
	end
	%b = ones(n);

	a_C = {};
	aa = [];
	for idx=1:nvar
	for jdx=1:nvar
		%a_C{idx,jdx} = dt*rand(n(1)*n(2),1);
		aa(idx,jdx,:) = dt*s*randn(n(1)*n(2),1);
	end
		%a_C{idx,idx} = 1 + a_C{idx,idx}; 
		aa(idx,idx,:) = 1 + aa(idx,idx,:);
	end

	% matrix setup
	if (n(1)<=128)
	I  = speye(prod(n));
	dx = L./n;
	[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,dx,2,{'circular','circular'},true);
	D2 = D2x+D2y;
	A = [];
	for idx=1:nvar
		A_ = [];
		for jdx=1:nvar
			if (jdx ~= idx)
				A_  = [A_,diag(sparse(squeeze(aa(idx,jdx,:))))];
			else
				A_  = [A_,diag(sparse(squeeze(aa(idx,jdx,:)))) - (dt*d(idx))*D2];
			end
		end
		A = [A;A_];
	end

	tic();
	x_ref = A \ b(:);
	toc()
	x_ref = reshape(x_ref,n(1),n(2),nvar);
	else
		x_ref = NaN(n(1),n(2),nvar);
	end

	%a0 = 1;

	if (0)
	mg = Multigrid();
	mg.init(a0,dt*d,L,n,1);
	
	tic()
	[x_mg,resn]=mg.mldivide(b,x0);
	toc()
	end


	mg_j = Multigrid_java();
	mg_j.init(aa,dt*d,L,n);
	mg_j.reltol = sqrt(eps);
	mg_j.nmaxiter = 100;	

	if (1)
	x_j = mg_j.mldivide(reshape(b,[],nvar)',reshape(x0,[],nvar)');
	x_j = reshape(x_j',n(1),n(2),nvar);
	relres = mg_j.relres
	iter = mg_j.iter
	
	if (n(1)<=128)
	rms(x_j-x_ref,'all')
	end
	%relres(end)
	%length(resn)
	end

	if (1)
	figure(1)
	clf();
	for idx=1:nvar
	subplot(nvar,3,1+(idx-1)*3);
	imagesc(x_ref(:,:,idx));
	colorbar
	drawnow
	
	subplot(nvar,3,2+(idx-1)*3);
	imagesc(x_j(:,:,idx));	
	colorbar
	drawnow	

	subplot(nvar,3,3+(idx-1)*3);
	imagesc(x_j(:,:,idx)-x_ref(:,:,idx));
	colorbar
	drawnow	
%	colorbar
	end
	end
%	tic(); x_g = gmres(A,b(:)); toc()

if (0)
	% bencmark
	rt = [];
	maxiter=(1:10)';

	for idx=1:length(maxiter)
	mg_j.nmaxiter = maxiter(idx);
	tic
	x_j = mg_j.mldivide(reshape(b,[],nvar)',reshape(x0,[],nvar)');
	x_j = reshape(x_j',n(1),n(2),nvar);
	rt(idx,1)=toc;
	end
	subplot(2,2,4)
	clf
	plot(maxiter,rt);
	A=[ones(size(maxiter)),maxiter];
	c=A \ rt
	hold on;
	plot(maxiter,A*c)

	m = 30;
	tic
	for idx=1:m
		mg_j.nmaxiter = 1;
		x_j = mg_j.mldivide(reshape(b,[],nvar)',reshape(x0,[],nvar)');
		x_j = reshape(x_j',n(1),n(2),nvar);
	end
	toc/m
	x_ = x_j;
	tic
	for idx=1:m
		mg_j.nmaxiter = 1;
		x_j = mg_j.cycle1(reshape(b,[],nvar)');
		%mg_j.cycle1(reshape(b,[],nvar)');
		%x_j = mg_j.
		x_j = reshape(x_j',n(1),n(2),nvar);
	end
	toc/m
	rms(x_-x_j,'all')

end
