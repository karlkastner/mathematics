% Wed 19 Jul 17:38:03 CEST 2023
function [z,Al,Ar] = step_advection_diffusion_trapezoidal(dt,dx,n,z,a,e)
	% TODO, throw error if stability conditions are not satisfied
	% TODO use krylov subspace
	D1 = (0.5*dt*a(1))*derivative_matrix_1_1d(n(1),dx(1),2,'circular','circular',true);
	D2 = (0.5*dt*e(1))*derivative_matrix_2_1d(n(1),dx(1),2,'circular','circular',true);
	I  = speye(n);
	if (length(n)>1)
		D1y = (0.5*dt*a(2))*derivative_matrix_1_1d(n(2),dx(2),2,'circular','circular',true);
		D2y = (0.5*dt*e(2))*derivative_matrix_2_1d(n(2),dx(2),2,'circular','circular',true);
		Iy  = speye(n(2));
		D1  = kron(D1,Iy) + kron(I,D1y);
		D2  = kron(D2,Iy) + kron(I,D2y);
		I   = speye(prod(n));
	end
	Al  = (I - D1 - D2);
	Ar  = (I + D1 + D2);
	z   = Al \ (Ar*z);
end

