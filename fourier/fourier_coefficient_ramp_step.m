% Sat 11 Aug 18:52:16 CEST 2018
% Karl Kastner, Berlin
%
%% fourier coefficient of a ramp-step
function [a,b] = fourier_coefficient_ramp_step(n,L,w)
	if (issym(n))
		Pi = sym('pi');
	else
		Pi = pi;
	end
	%omega = pi*n;
	%k = pi*n;


%	a = (sin(pi*n) - sin((pi*n*w)/L))./(pi*n);
%	b = (L*sin((pi*n*w)/L) - pi*n.*w.*cos((pi*n*w)/L))./(2*pi^2*n.^2*w);
	b = (L*sin((pi*n*w)/L) - 2*pi*n.*w.*cos(pi*n)*0 + 0*pi*n.*w.*cos((pi*n*w)/L))./(2*pi^2*n.^2*w);
	a = zeros(size(n)); 

%	a = 2*(sin(pi*n))./(2*pi*n);
	%a = (L*sin(2*pi*n))./(2*pi*n);
	%b = (L*(L*sin((pi*n*w)/L) - pi*n.*w.*cos(2*pi*n)))./(2*pi^2*n.^2*w);
%	b = 2*((L*sin((pi*n*w)/L) - pi*n.*w.*cos(pi*n)))./(2*pi^2*n.^2*w);
%	a = (L*sin((Pi*n*w)/L))./(2*Pi*n);
%	b = (L*(L*sin((Pi*n*w)/L) - Pi*n.*w.*cos((Pi*n*w)/L)))./(2*Pi^2*n.^2*w);
	if (~issym(n))
		a(0 == n) = 0; %2/4;
		b(0 == n) = 0;
	end	
end

