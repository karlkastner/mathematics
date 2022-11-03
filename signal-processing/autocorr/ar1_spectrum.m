% Fri  2 Feb 16:20:39 CET 2018
%% spectrum of the ar1 process
function [H] = ar1_spectrum(rho,n)
	if (1) % continuous
		L = n;
		x = (0:L-1)'/L;
		%H = (1+rho^2-2*rho)./(1+rho^2-2*rho*cos(2*pi*x));
		H = (1+rho^2-2*rho)./(1+rho^2-2*rho*cos(2*pi*x));
	else % discrete
		% TODO somehow the change from rho to a does not work properly
		x = (1:n)'/n;
		a= 1-rho;
		Ts= 1;
		z = exp(1i*2*pi*x);
		% Laplace: Hp(s) := a/(s + a)
		% y_i = (1-exp(-aT)) x_i +  exp(-aT) y_i-1
		% => (1-rho)rho = exp(-aT)
		a = -1/Ts*log(rho*(1-rho));
		H(:,2) = abs( (1-exp(-a*Ts))./(1-exp(-a*Ts)./z) );
	end
end

