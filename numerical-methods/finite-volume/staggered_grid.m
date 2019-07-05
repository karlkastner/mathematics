% 2013-04-20 16:52:50.000000000 +0200
%% staggered grid approximation to the SWE

function [h_ m_] = staggered_grid(dt,dx,h,m)
	m_l =
	m_r =

	h_c = 

	B = 	
	
	m_c = 0.5*(m(1:end-1) + m(2:end));
	h_c = 0.5*(h(1:end-1) + h(2:end));

	[A_l B] = feval(qfun, [], m(1:end-1));
	[A_r B] = feval(qfun, [], m(2:end));

	[A B_l] = feval(qfun, h(1:end-1), 0.5*(m(1:end-1) + m(2:end-1)));
	[A B_r] = feval(qfun,   h(2:end), 0.5*(m(2:end-1) + m(3:end));
%	[A B_l] = feval(qfun, h(1:end-1), sqrt(2)*m(1:end-1)));
%	[A B_r] = feval(qfun,   h(2:end), sqrt(2)*m(2:end-1));

	h = h + 0.5*dt/dx*(A_r - A_l);
	m = m + 0.5*dt/dx*(B_r - B_l);
end

