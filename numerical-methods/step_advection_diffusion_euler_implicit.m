% Wed 19 Jul 17:38:03 CEST 2023
function [z,A] = step_advection_diffusion_implicit_euler(dt,dx,n,z,a,e)
	% TODO use krylov
	D1x = derivative_matrix_1_1d(n,dx,2,'circular','circular',true);
	D2x = derivative_matrix_2_1d(n,dx,2,'circular','circular',true);
	Ix  = speye(n);
	A = (Ix - (dt*a)*D1x - (dt*e)*D2x);
	z = A \ z;
end

