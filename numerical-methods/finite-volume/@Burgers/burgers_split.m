% Fri Apr 22 22:23:59 MSD 2011
% Karl Kästner
% 
%% viscous Burgers' equation,
%% mixed analytic and numerical derivative in frequency space
%% by splitting sheme
%% u_t = -(0.5*u^2)_x + c*u_xx
function u = burgers_split(u, dt, c, Df1, Df2)
	% linear part - analytic
	u_fft = fft(u);
	u_fft = exp(c*dt*Df2).*u_fft;
	u = ifft(u_fft);
	% non-linear part - explicit
	u = u - u.*ifft(Df1.*u);
end % burgers_split

