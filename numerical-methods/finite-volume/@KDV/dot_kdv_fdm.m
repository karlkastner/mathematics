% Sat Apr 23 22:02:25 MSD 2011
% Karl Kästner
%
%% korteweg de vries equation
%% u_t + (0.5*u^2)_x = c*u_xxx
function u_dot = kdv(t, u, c, D1, D3)
	u_dot = -D1*0.5*u.^2 + c*D3*u;
	%u_hat = fft(u);
	%F = 6*u.*ifft(k.*u_hat) + ifft(kcub.*u_hat);
end % kdvFft

