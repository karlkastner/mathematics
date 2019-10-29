% Wed Jun 27 16:10:03 MSK 2012
% Karl Kästner, Berlin

function A = poisson(n)
	switch(length(n))
		case {1}
			A = (n(1)+1)^2*spdiags(ones(n(1),1)*[1 -2 1],-1:1,n(1),n(1));
		case {2}
			A1 = (n(1)+1)^2*spdiags(ones(n(1),1)*[1 -2 1],-1:1,n(1),n(1));
			A2 = (n(2)+1)^2*spdiags(ones(n(2),1)*[1 -2 1],-1:1,n(2),n(2));
			I1=speye(n(1));
			I2=speye(n(2));
			A = kron(A1,I2) + kron(I1,A2);
		case {3}
			A1 = (n(1)+1)^2*spdiags(ones(n(1),1)*[1 -2 1],-1:1,n(1),n(1));
			A2 = (n(2)+1)^2*spdiags(ones(n(2),1)*[1 -2 1],-1:1,n(2),n(2));
			A3 = (n(3)+1)^2*spdiags(ones(n(3),1)*[1 -2 1],-1:1,n(3),n(3));
			I1=speye(n(1));
			I2=speye(n(2));
			I3=speye(n(3));
			A = kron(kron(A1,I2),I3) + kron(kron(I1,A2),I3) + kron(kron(I1,I2),A3);
		case {4}
			A1 = (n(1)+1)^2*spdiags(ones(n(1),1)*[1 -2 1],-1:1,n(1),n(1));
			A2 = (n(2)+1)^2*spdiags(ones(n(2),1)*[1 -2 1],-1:1,n(2),n(2));
			A3 = (n(3)+1)^2*spdiags(ones(n(3),1)*[1 -2 1],-1:1,n(3),n(3));
			A4 = (n(4)+1)^2*spdiags(ones(n(4),1)*[1 -2 1],-1:1,n(4),n(4));
			I1=speye(n(1));
			I2=speye(n(2));
			I3=speye(n(3));
			I4=speye(n(4));
			A = kron(kron(kron(A1,I2),I3),I4) + kron(kron(kron(I1,A2),I3),I4) + kron(kron(kron(I1,I2),A3),I4)+kron(kron(kron(I1,I2),I3),A4);
	end
end % function poisson

