% Sat 11 Aug 18:52:16 CEST 2018
% Karl Kastner, Berlin
% 
%% fourier series coefficient of a ramp pules
function [a,b] = fourier_coefficient_ramp_pulse(n,L,w)
	if (issym(n))
		Pi = sym('pi');
	else
		Pi = pi;
	end
	%omega = pi*n;
	%k = pi*n;
	a = (L*sin((Pi*n*w)/L))./(2*Pi*n);
	b = (L*(L*sin((Pi*n*w)/L) - Pi*n.*w.*cos((Pi*n*w)/L)))./(2*Pi^2*n.^2*w);
	if (~issym(n))
		a(0 == n) = w/4;
		b(0 == n) = 0;
	end	
end

