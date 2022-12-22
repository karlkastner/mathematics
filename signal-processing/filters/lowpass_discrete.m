% Fri  4 Aug 15:46:21 CEST 2017
%
%% design coefficients of a low pass filter with specified cut of frequency
%% and sampling period
%% alalogue low pass with pole at s=-omega_c=1/tau=1/RC
%% Ha = tau/(tau + s) = 1/(1 + omega_c*s)
%
function [b a] = digital_low_pass_coeficients(fc,dt)
	% angular frequency
	omega_c = 2*pi*fc;

	% warp frequency
	omega_c = (2/dt)*tan(omega_c*dt/2);
	
%	omega_c = 1/omega_c;

	% low-pass analogue to digital
	% c.f. MIT Lec. 19
	b = (1-exp(-omega_c*dt));
	a = [1,-exp(-omega_c*dt)];
end

