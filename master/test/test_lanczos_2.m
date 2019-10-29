function test_lanczos()

%f = {@setup_1D_, @setup_2D_}
f = {@setup_1D_, @setup_2D_, @setup_3D_}
%f = {@setup_1D, @setup_2D, @setup_3D}
%f = {@setup_2D}

for fdx=1:length(f)
	% setup the matrix
	A = feval(f{fdx});
	[Q T q_mp1 beta_mp1 alpha beta] = lanczos(A,length(A));

	E_ = sort(eig(A));

%	sort(E_)
%	pause

	tol = sqrt(eps);

	K(1,1)=0;
	for mdx=2:length(T)
		[S Theta] = eigs(T,mdx);
		Theta = diag(Theta);
		% select converged eigenvalues
		T2 = [];
		for idx=1:length(Theta)
			if (beta(mdx+1)*abs(S(end,idx))/abs(Theta(idx)) < tol)
				T2(end+1,1) = Theta(idx);
			end
		end
		% sort eigenvalues
		T2 = sort(T2);
		% remove duplicates
		T3 = [];
		if (length(T2) > 1)
			T3(1,1) = T2(1,1);
		end
		for idx=2:length(T2)
			if (abs(T2(idx,1) - T2(idx-1,1))/abs(T2(idx-1,1)) > tol)
				T3(end+1,1) = T2(idx,1);
			end
		end
		K(mdx,1) = length(T3);
		Err(mdx,1) = norm(Q(:,1:mdx)'*Q(:,1:mdx));
		if (length(T3) > 1)
			Err_(mdx,1) = norm(T3 - E_(1:length(T3)));
		end
		[mdx length(T2) length(T3) Err(mdx)]
	end

	figure(fdx)	
	subplot(2,2,1)
	plot(K,'k','Linewidth',2); % Err Err_])
	grid on                                  
	xlabel('iteration')
	ylabel('number of converged eigenvalues')
	if (1==fdx)
		%title('Convergence of the Eigenvalues of the 1D-Laplacian by Lanczos Iteration')
		title('1D-Laplacian') %Convergence of the Eigenvalues of the 1D-Laplacian by Lanczos Iteration')
		xlim([0 255])
		ylim([0 260])
		print -deps ../img/lanczos-iteration-1d.eps
	else
		xlim([0 100])
		%title('Convergence of the Eigenvalues of the 2D-Laplacian by Lanczos Iteration')
		title('2D-Laplacian') % by Lanczos Iteration')
		print -deps ../img/lanczos-iteration-2d.eps
	end
end % fdx

end % test lanczos 2

function A = setup_1D_()
	n = 100;
	A = harmonic_oscillator(n,1);
end

function A = setup_2D_()
	n = 10;
	A = harmonic_oscillator(n, 2, 1);
end

function A = setup_3D_()
	n = 5;
	A = harmonic_oscillator(n, 3, 1, 1);
end

function A = setup_1D()
%	n = 250;
	n = 100;
	A = spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
end

function A = setup_2D()
	%n1 = 10; %20;
	%n2 = 25; %50;
	n1 = 10;
	n2 = 10;

	% setup the matrix
	L1 = spdiags(ones(n1,1)*[1 -2 1],-1:1,n1,n1);
	L2 = spdiags(ones(n2,1)*[1 -2 1],-1:1,n2,n2);
	L2 = 0.9999*spdiags(ones(n2,1)*[1 -2 1],-1:1,n2,n2);
	I1 = speye(n1);
	I2 = speye(n2);
	A = kron(L1,I2) + kron(I1,L2);
end

function A = setup_3D()
	%n1 = 10; %20;
	%n2 = 25; %50;
	n1 = 5;
	n2 = 5;
	n3 = 5;

	% setup the matrix
	L1 = 1*spdiags(ones(n1,1)*[1 -2 1],-1:1,n1,n1);
	L2 = 0.9999*spdiags(ones(n2,1)*[1 -2 1],-1:1,n2,n2);
	L3 = 0.9999^2*spdiags(ones(n3,1)*[1 -2 1],-1:1,n3,n3);
	I1 = speye(n1);
	I2 = speye(n2);
	I3 = speye(n3);
	A = kron(kron(L1,I2),I3) + kron(kron(I1,L2),I3) + kron(kron(I1,I2),L3);
end

