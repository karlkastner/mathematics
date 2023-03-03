% 2023-01-20 09:36:10.846094957 +0100
syms fx sx f0 c
pi_ = sym('pi');
	z = 1i*(2*f0*fx*pi_*sx^2) + (pi_^2*f0^2*sx^4 + f0^2 - fx.^2);    
	solve(exp(1i*pi*c) - z, fx)
