% Sun Apr 24 14:28:39 MSD 2011
% Karl Kästner
% 
%% viscous burgers' equation
%% u_t = -d/dx (1/2*u^2) + c d^2/dx^2 u_xx
% TODO why fft?
function u_dot = burgers(t, u, c, D1, D2)
	u_dot = -D1*fft(0.5*u.^2) + c*D2*u;
end % burger

