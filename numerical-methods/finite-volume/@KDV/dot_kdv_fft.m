% Fri Apr 22 22:23:28 MSD 2011
% Karl Kästner
%
%% korteweg de vries equation
%% compute derivatives in frequency space
%% u_t + (0.5*u^2)_x = c*u_xxx
% TODO work with real inputs
function u_dot = dot_kdv_fft(t, u_fft, c, D1, D3)
	u     = ifft(u_fft);
	u_dot = -6*D1.*fft(0.5*u.^2) + c*D3.*u_fft; % 0.5 ?
end % dot_kdv_fft

