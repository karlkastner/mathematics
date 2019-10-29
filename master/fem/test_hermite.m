% Thu Feb 14 02:47:37 MSK 2013
% Karl Kästner, Berlin

%lagrange vs hermite element matrices (1D):
%y=2.^(-15:5)';
%	for idx=1:length(y); x=y(idx);
%	A = [	1 0 0 0;
%		1 x x^2 x^3;
%		0 1 0 0;
%		0 1 2*x 3*x^2];
%	C(idx,1)=cond(A);
%	A = [1 0  0 0;
%		1 1/3*x 1/9*x^2 1/27*x^3;
%		1 2/3*x 4/9*x^2 8/27*x^3;
%		1 x x^2 x^3];
%	C(idx,2)=cond(A); end; loglog(y,C), legend('hermite','lagrange'); grid on

function test_hermite()
	N = 2.^(0:7)+2;
	addpath('/home/pia/Documents/master/thesis/src/fem/int/');

	v_func = {'lagrange'; 'hermite'}

	for vdx=1:length(v_func);

	for ndx=1:length(N)
		n = N(ndx);
		
		% points
		P = linspace(0,1,n)';
		% elements
		T = [ (1:n-1)' (2:n)' ];
	
		% higher order
		if (1 == vdx)
			h = 1/(n+1);
			P = [P; P+1/3*h; P+2/3*h];
		else
			P = [P; P];
		end
		T = [T T+n];
	
		[w b] = feval('int_1d_gauss_6');
	
		% assemble mass and stiffness matrix
		A = assemble_dphi_dphi(P,T,v_func{vdx},b,w);
		B = assemble_phi_phi(P,T,v_func{vdx},b,w);
	
		% apply homogenous dirichlet boundary conditions
		A(:,n) = []; A(n,:) = [];
		B(:,n) = []; B(n,:) = [];
		A(:,1) = []; A(1,:) = [];
		B(:,1) = []; B(1,:) = [];
		
		% compute the smallest eigenvalue
		e = eigs(A,B,1,'SM');
		E(ndx,vdx) = e;
	end % for ndx
	end % for vdx
	abserr = E-pi^2
	(diff(abserr')./sum(abserr'))'
	relerr = abserr/pi^2;
	loglog(N,abs(relerr))
	grid on
	legend(v_func)
end % test hermite
	
	

function A = assemble_phi_phi(P,T,v_func,b,w)

	m = 0;
	buf = [ 0 0 0];
	for tdx=1:size(T,1)
		% triangle points
		A = [ones(size(T,2),1) P(T(tdx,:))];

		% quadrature points
		q = b*A(1:2,1:2);

		% vandermonde matrix
		A = feval(v_func,A(:,2));

		% test functions
		C = inv(A);

		% quadrature point vandermonde matrix
		Vq = [q q(:,2).^2 q(:,2).^3];

		% evaluate the test functions at the quadrature points
		phi = Vq*C;
	
		% matrix entries
		for idx=1:4
			for jdx=1:4
				I = sum(w.*phi(:,idx).*phi(:,jdx));
				m = m+1;
				buf(m,:) = [T(tdx,jdx) T(tdx,idx) I];
			end % for jdx
		end % for idx
	end % for tdx
	A = sparse(buf(:,1), buf(:,2), buf(:,3));
end % function assemple_phi_phi

function A = assemble_dphi_dphi(P,T,v_func,b,w)

	m = 0;
	buf = [ 0 0 0];
	for tdx=1:size(T,1)
		% triangle points
		A = [ones(size(T,2),1) P(T(tdx,:))];

		% quadrature points
		q = b*A(1:2,1:2);

		% vandermonde matrix
		A = feval(v_func, A(:,2));

		% test functions
		C = inv(A);

		% test function derivatives
		dC = [ C(2,:); 2*C(3,:); 3*C(4,:) ];

		% quadrature point vandermonde matrix
		 %Vq = feval(v_func,q);
		 %Vq = Vq(:,1:3);
		Vq = [q q(:,2).^2];

		% evaluate derivatives
		dPhi = Vq*dC;

		% assemble matrix
		for idx=1:4
			for jdx=1:4
				I = sum(w.*dPhi(:,idx).*dPhi(:,jdx));
				m = m+1;
				buf(m,:) = [T(tdx,idx) T(tdx,jdx) I];
			end % for jdx
		end % for idx
	end % for tdx
	A = sparse(buf(:,1), buf(:,2), buf(:,3));
end % function assemble_dphi_dphi

function A = lagrange(x)
	A = [	 1, x(1) x(1)^2 x(1)^3;
		 1, x(2) x(2)^2 x(2)^3;
		 1, x(3) x(3)^2 x(3)^3;
		 1, x(4) x(4)^2 x(4)^3];
end

function A = hermite(x)
		A = [	1   x(1)   x(1)^2     x(1)^3;
			1 x(2) x(2)^2   x(2)^3;
			0      1   2*x(1)   3*x(1)^2;
			0      1 2*x(2) 3*x(2)^2];
end

