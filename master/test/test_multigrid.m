
%n=10; A = poisson([n n n]); E = eig(full(A)); A=A+1/12*A^2/(n+1)^2; H=eig(full(A)); [A B] = poisson_fem([n n n]); F = eig(full(A),full(B)); G = -pi^2*(length(F):-1:1)'.^2; plot([E F 0.5*(E+F) H]); %semilogy(abs([E-G 0.5*(E+F)-G F-G]))  
%n=100; A = poisson([n]); E = eig(full(A)); A=A+1/12*A^2/(n+1)^2; H=eig(full(A)); [A B] = poisson_fem([n]); F = eig(full(A),full(B)); G = -pi^2*(length(F):-1:1)'.^2; plot([E F 0.5*(E+F) G H]); %semilogy(abs([E-G 0.5*(E+F)-G F-G])) 
%n=100000; A = poisson(n); tic(); eigs(A,[],1,'SM'); toc(), [A B] = poisson_fem(n); tic(); eigs(A,B,1,'SM'); toc()
%n=100; A = poisson([n n]); tic(); eigs(A,[],10,'SM'); toc(), [A B] = poisson_fem([n n]); tic(); eigs(A,B,10,'SM'); toc()

function test_multigrid()
	N = unique(ceil(2.^(0:0.5:5.5)));
	for idx=1:length(N)
		n = N(idx)
		A = poisson([n n]);
		b = rand(size(A,1),1);
		x = zeros(size(A,1),1);
		tol = 1e-7;
	        [x flag nrvec r] = jacobi(A,b,x,tol);
		M(idx) = length(nrvec)
		%[flag length(nrvec)]
	end
	loglog(N.^2,M)
end


function [x flag nrvec r] = v_cycle(Ac, Mc, b, x, tol)
	
	x = v_step(Ac, Mc, b_, x_, idx);
		
	end
end
	
