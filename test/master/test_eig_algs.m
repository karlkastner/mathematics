function test_eig_algs()
	d = 3;
	n(1) = 68;
	n(1) = 20;
	n(2) = n(1)/2;
	n(3) = n(1)/4;
	I1 = speye(n(1));
	A1 = (n(1)+1)^2*spdiags(ones(n(1),1)*[1 -2 1],-1:1,n(1),n(1));
	I2 = speye(n(2));
	A2 = (n(2)+1)^2*spdiags(ones(n(2),1)*[1 -2 1],-1:1,n(2),n(2));
	I3 = speye(n(3));
	A3 = (n(3)+1)^2*spdiags(ones(n(3),1)*[1 -2 1],-1:1,n(3),n(3));
	switch(d)
		case{1}
			A = A1;
			x0 = rand(n(1),1);
		case {2}
			A = kron(A1,I2) + kron(I1,A2);
			x0 = rand(n(1)*n(2),1);
		case {3}
			A = kron(kron(A1,I2),I3) + kron(kron(I1,A2),I3) + kron(kron(I1,I2),A3);
			n_ = n(1)*n(2)*n(3);
			x0 = rand(n(1)*n(2)*n(3),1);
	end

	tic()
	eigs(A,[],10);
	toc()

	Z0 = rand(n_,10);
	tic()
	E = simultaneous_iteration(A,Z0);
	toc()

	tic()
	eigs(A,[],10,'SM');
	toc()

	Z0 = rand(n_,10);
	tic()
	E = simultaneous_inverse_iteration(A,Z0);
	toc()

	
	return

%	x0 = sin(pi/n_*(1:n_))';
%	x1 = sin(2*pi/n_*(1:n_))';
%	[x flag_ relres iter1] = minres(A,x0,[],n(1)*n(2)*n(3),[],[],[]);
%	iter1
%	%rand(n(1)*n(2)*n(3),1);
%	x1_ = x1 - (x1'*x0)/(x0'*x0)*x0;
%	norm(x0)
%	norm(x1)
%	norm(x1_)
%	[x flag_ relres iter1] = minres(A,x1_,[],n(1)*n(2)*n(3),[],[],[]);
%	iter1
%	x1_'*x0
%	pause 
	

	%eigs(A,[],6,'SM')

%{
	e_ = eigs(A,[],1);
	[e x1 L] = power_iteration(A,x0,0);
	disp(sprintf('%f %f %d',e,e-e_,length(L)))
	semilogy(abs(L-e_)); hold on

	% obviously, the power iteration always converges to the largest eigenvalue
	% even if if the start vector is orthogonal to it (just delays convergence)
	x0 = x0 - (x1'*x0)*x1;
	[e x2 L] = power_iteration(A,x0,0);
	disp(sprintf('%f %f %d',e,e-e_,length(L)))
	semilogy(abs(L-e_),'r'); hold off
%}
	tic()
	e_ = eigs(A,[],1,'SM');
	toc()
	tic()
	[e x1 L] = inverse_iteration(A,x0,0);
	disp(sprintf('%f %f %f %d %f',e_,e,e-e_,length(L), toc()))
	semilogy(abs(L-e_)); hold on
	tic()
	[e x1 L] = inverse_iteration(A,x0,1);
	disp(sprintf('%f %f %f %d %f',e_,e,e-e_,length(L), toc()))
	tic()
	[e x1 L] = inverse_iteration(A,x0,2);
	disp(sprintf('%f %f %f %d %f',e_,e,e-e_,length(L), toc()))
	tic()
	[e x1 L] = inverse_iteration(A,x0,3);
	disp(sprintf('%f %f %f %d %f',e_,e,e-e_,length(L), toc()))
	eigs(A,[],6,'SM')

%	x0 = x0 - (x1'*x0)*x1;
%	e_ = eigs(A,[],1,'SM');
%	[e x2 L] = inverse_iteration(A,x0,1);
%	disp(sprintf('%f %f %f %d',e_,e,e-e_,length(L)))
%	semilogy(abs(L-e_),'r'); hold off

	
end


function E = simultaneous_iteration(A,Z)
	tol = sqrt(eps);
	[Q R] = qr(Z,0);
	edx=0;
	sigma=R(edx+1,edx+1);
	E_ = eigs(A,10);
	idx=1;
	while (1)
		Z = A*Q;
		[Q R] = qr(Z,0);
		sigma_old = sigma;
		sigma=R(edx+1,edx+1);
		% check for convergence
		idx=idx+1;
		if (abs(sigma_old - sigma) < tol)
			% save the converged eigenvalue
			edx=edx+1;
			E(edx,1) = R(edx,edx);
			% check if all eigenvalues converged
			%if (1==size(Q,2))
			if (edx==size(Q,2))
				return;
			end
			% deflate - no deflation, otherwise only convergence to maximum eigenvalue!
			%Q = Q(:,2:end);
		end
	end % while 1
end % simultaneous iteration

function E = simultaneous_inverse_iteration(A,Z)
	tol = sqrt(eps);
	[Q R] = qr(Z,0);
	I = speye(size(A));
	edx=0;
	sigma=1/R(edx+1,edx+1);
	E_ = eigs(A,10);
	idx=1;
        [L_,U_,pp_,qq_,D_] = lu(A);
	while (1)
		%Z = (A - sigma*I) \ Q;
		Z = qq_*(U_ \ (L_ \ (pp_*(D_\Q))));
		[Q R] = qr(Z,0);
		sigma_old = sigma;
		sigma=1/R(edx+1,edx+1); % + sigma;
		% check for convergence
		idx=idx+1;
		if (abs(sigma_old - sigma) < tol)
			% save the converged eigenvalue
			edx=edx+1;
			E(edx,1) = sigma;
			% check if all eigenvalues converged
			%if (1==size(Q,2))
			if (edx==size(Q,2))
				return;
			end
			% deflate - no deflation, otherwise only convergence to maximum eigenvalue!
			%Q = Q(:,2:end);
		end
	end % while 1
end % simultaneous iteration

function [lambda x L] = inverse_iteration(A,x,flag)
	tol = 1e-7;
	maxiter = 10*size(A,1);
	I = speye(size(A));
	x = 1/sqrt(x'*x)*x;
	lambda = 0;
	idx=0;
	if (2 == flag)
                [L_,U_,pp_,qq_,D_] = lu(A);
	end
               %[L,U,pp,qq,D] = lu(A);
	while (1)
		x_old = x;
		switch(flag)
			case {0} % shift and exact
				sigma = lambda;
				x = (A - sigma*I) \ x;
			case {1} % shift and not exact
				x0=x;
				sigma = lambda;
				[x flag_ relres iter1] = minres(A-sigma*I,x,tol,maxiter,[],[],x0);
			case {2} % no shift and exact
				sigma=0;
                		x = qq_*(U_ \ (L_ \ (pp_*(D_\x))));
			case {3} % no shift and no exact
				sigma=0;
				x0 = x;
				[x flag_ relres iter1] = minres(A,x,tol,maxiter,[],[],x0);
		end
		lambda_old = lambda;
		lambda = 1/(x'*x_old) + sigma;
		x = 1/sqrt(x'*x)*x;
		idx = idx+1;
		L(idx) = lambda;
		if ( abs(lambda - lambda_old) < tol || idx > maxiter)
			return;
		end
	end
end

