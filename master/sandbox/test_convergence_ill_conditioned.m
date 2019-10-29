% Tue Mar 12 02:08:47 MSK 2013
% Karl Kästner, Berlin

% impact on refining the uniform mesh on convergence of the discrete laplacian
% first diemension is refined twice as fast and the limit matrix is ill-conditioned

function test_convergence_ill_conditioned()

	% finding : both series converge, in 2d:
	%	    round off error of ill conditioned mesh is just twice as large and therefore irrelevant until maximum accuracy reached (err = eps*cond(A))
	%           ill conditioned mesh is even twice as accurate (as effectively only one dimension contributes to the error)
	%		=> round of determined by smallest element, discritisation error determined by largest element

	n(1) = 1;
	n(2) = 1;

	d=2;

	opts.tol = 1e-12;
	order = 10;

	tic()
	for idx=1:5 % 6

	A1 = (1+n(1))^2*spdiags(ones(n(1),1)*[1 -2 1],-1:1,n(1),n(1));
	I1 = speye(n(1));
	A2 = (1+n(2))^2*spdiags(ones(n(2),1)*[1 -2 1],-1:1,n(2),n(2));
	I2 = speye(n(2));
	if (2 == d)
		A3 = 1; I3 = 1; n(3)=1;
	else
		A3 = A1; I3 = I1; n(3) = n(1);
	end

	A1_ = A1;
	A2_ = A2;
	for jdx=2:order
		jdx
		jdx_=jdx-1
		c = nchoosek(2*jdx_+1,jdx_)*(jdx_+1)^2
		A1 = A1 - 1/c*(-1)^jdx*1/(n(1)+1)^(2*jdx-2)*A1_^jdx;
		A2 = A2 - 1/c*(-1)^jdx*1/(n(2)+1)^(2*jdx-2)*A2_^jdx;
	end

	A = kron(kron(A1,I2),I3) + kron(kron(I1,A2),I3) + (d>2)*kron(kron(I1,I2),A3);
	I = speye(n(1)*n(2)*n(3));

	size(A)
	e(idx,1) = eigs(A,[],1,'SM',opts)
	f(idx,1) = eigs(A,[],1,'LM',opts)

	A = kron(kron(A1,I1),I3) + kron(kron(I1,A1),I3) + (d>2)*kron(kron(I1,I1),A3);
	I = speye(n(1)^2*n(3));
	size(A)
	e(idx,2) = eigs(A,[],1,'SM',opts)
	f(idx,2) = eigs(A,[],1,'LM',opts)
	
	N(idx,:) = n;
	n(1)=n(1)*2;
	n(2)=n(2)*4;
	
	end
	toc()
	e
	err = e + d*pi^2
	e_diff = abs(diff(e')');
	frr = e./f
	N
	semilogy([abs(err) frr (N(:,1)./N(:,2)).^2],'.-');

end


