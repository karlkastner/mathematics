% Wed Apr 27 00:22:23 MSD 2011
% Karl Kästner
%
%% korteweg de vries equation in frequency space,
%% derivative treated by splitting scheme
function u = kdv_split(u, dt, Df1, Df3)
	% linear part - analytic
	u_hat = fft(u);
	u_hat = exp(dt*Df3).*u_hat;
	% linear part - implicit
	%u_hat = (1 - dt*Df3).^-1.*u_hat;
	% non-linear part - explicit
	u = ifft(u_hat);
	u = u + dt.*6*u.*ifft(Df1.*u_hat);
	% u = exp(dt*(6*Df1.*u + D3).*u;
end % kdv_split

