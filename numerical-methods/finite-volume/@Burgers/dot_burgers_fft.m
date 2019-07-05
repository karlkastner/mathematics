% Fri Apr 22 22:23:59 MSD 2011
% Karl Kästner
% 
%% viscous Burgers' equation in frequency space
%% u_t + (0.5*u^2)_x = c*u_xx
function u_dot = dot_burgers_fft(t, u_fft, c, D1, D2)
	u     = ifft(u_fft);
	u_dot = -D1.*fft(0.5*u.^2) + c*D2.*u_fft;
end % burgerFft

