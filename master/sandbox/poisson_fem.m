% Wed Jul 11 13:46:02 MSK 2012
% Karl Kästner, Berlin

% discretise the laplacian operator with pice wise linear lemenent
% on a homogenous grid with d-dimensional hypercubes

% n=100; A = poisson(n); E = eig(full(A)); [A B] = poisson_fem(n); F = eig(full(A),full(B)); G = -pi^2*(length(F):-1:1)'.^2; plot([E F G]); semilogy(abs([E-G F-G])
%n=10; A = poisson([n n n]); E = eig(full(A)); A=A+1/12*A^2/(n+1)^2; H=eig(full(A)); [A B] = poisson_fem([n n n]); F = eig(full(A),full(B)); G = -pi^2*(length(F):-1:1)'.^2; plot([E F 0.5*(E+F) H]); %semilogy(abs([E-G 0.5*(E+F)-G F-G]))  
%n=100; A = poisson([n]); E = eig(full(A)); A=A+1/12*A^2/(n+1)^2; H=eig(full(A)); [A B] = poisson_fem([n]); F = eig(full(A),full(B)); G = -pi^2*(length(F):-1:1)'.^2; plot([E F 0.5*(E+F) G H]); %semilogy(abs([E-G 0.5*(E+F)-G F-G])) 
%n=100000; A = poisson(n); tic(); eigs(A,[],1,'SM'); toc(), [A B] = poisson_fem(n); tic(); eigs(A,B,1,'SM'); toc()
%n=100; A = poisson([n n]); tic(); eigs(A,[],10,'SM'); toc(), [A B] = poisson_fem([n n]); tic(); eigs(A,B,10,'SM'); toc()

% this is obviously not true in higher dimensions, try the following:
%	n = 100;
%	int = @int_1d_trapezoidal;
%	int = @int_1d_gauss_6;
%	[P T BC] = mesh_1d_uniform(n,1);
%	[P T BC] = promote_1d_2_3(P, T, BC);
%	A = assemble_1d_dphi_dphi(P,T,int,[]);
%	B = assemble_1d_phi_phi(P,T,int,[]);
%	[A B] = boundary_1d(A,B,BC);
%	eigs(A,B,6,'SM')
%	I = speye(n);
%	A=kron(A,I)+kron(I,A);
%	B=kron(B,I)+kron(I,B);
%	eigs(A,B,6,'SM')

function [A B] = poisson_fem(n)
	d = length(n);
	switch (d)
		case {1}
			A1 = (n(1)+1)*spdiags(ones(n(1),1)*[1 -2 1],-1:1,n(1),n(1));
			B1 = 1/(6*(n(1)+1))*spdiags(ones(n(1),1)*[ 1 4  1],-1:1,n(1),n(1));
			A = A1;
			B = B1;	
		case {2}
			A1 = (n(1)+1)*spdiags(ones(n(1),1)*[1 -2 1],-1:1,n(1),n(1));
			B1 = 1/(6*(n(1)+1))*spdiags(ones(n(1),1)*[ 1 4  1],-1:1,n(1),n(1));
			I1 = speye(n(1));
			A2 = (n(2)+1)*spdiags(ones(n(2),1)*[1 -2 1],-1:1,n(2),n(2));
			B2 = 1/(6*(n(2)+1))*spdiags(ones(n(2),1)*[ 1 4  1],-1:1,n(2),n(2));
			I2 = speye(n(2));
			A = kron(A1,I2) + kron(I1,A2);
			B = kron(B1,I2) + kron(I1,B2);
		case {3}
			A1 = (n(1)+1)*spdiags(ones(n(1),1)*[1 -2 1],-1:1,n(1),n(1));
			B1 = 1/(6*(n(1)+1))*spdiags(ones(n(1),1)*[ 1 4  1],-1:1,n(1),n(1));
			I1 = speye(n(1));
			A2 = (n(2)+1)*spdiags(ones(n(2),1)*[1 -2 1],-1:1,n(2),n(2));
			B2 = 1/(6*(n(2)+1))*spdiags(ones(n(2),1)*[ 1 4  1],-1:1,n(2),n(2));
			I2 = speye(n(2));
			A3 = (n(2)+1)*spdiags(ones(n(2),1)*[1 -2 1],-1:1,n(2),n(2));
			B3 = 1/(6*(n(2)+1))*spdiags(ones(n(2),1)*[ 1 4  1],-1:1,n(2),n(2));
			I3 = speye(n(2));
			A = kron(kron(A1,I2),I3) + kron(kron(I1,A2),I3) + kron(kron(I1,I2),A3);
			B = kron(kron(B1,I2),I3) + kron(kron(I1,B2),I3) + kron(kron(I1,I2),B3);
	end
end

