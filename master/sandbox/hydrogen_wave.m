% Tue Sep  6 17:25:47 MSD 2011
% Karl Kästner, Berlin

function psi = hydrogen_wave(r, theta, phi, n, l, m)

	% general constants

	% elementary charge
	e_electron = 1.60217656535e-19;
	% mass of resting electron [kg]
	m_electron = 9.1093829140-31;
	% plack constant
	h_planck = 6.6260695729e-34;
	% reduced plack constant
	h_bar = 0.5*h_planck/pi;
	% vacuum light speed [m/s]
	c_light = 2.99792458e8;
	% vacuum permittivity
	eps_0 = 8.854187817620e-12;
	% fine structure constant
	a0_fine = e_electron^2/(2*eps_0*h_planck*c_light);

	%%% Hydrogen atom
	
	% Bohr radius
	a0 = 4*pi*eps_0*h_bar^2/(m_electron*e_electron^2);
	% ?
	rho = 2*r/(n*a0);

	% laguerre polynomial
	L = laguerre(n-l-1, 2*l+1, rho);
	% spherical harmonics
	Y = shperical_harmonics(l,m, theta, phi);

	% exp or e ?
	psi = sqrt( (2/(n*a0))^3 * fak( n - l - 1)/(2*n*fak(n + l)^3) ) ...
		* exp(-0.5*rho)*pho^l * L * Y;
end

function Y = spherical_harmonics(l, m, theta, phi)
	N = norm_factor(l,m) 
	P = legendre(l,m,cos(theta));
	Y = 1.0/sqrt(2*pi)*N*P*exp(i*m*phi);
end

function N = norm_factor(l, m)
	N = sqrt(0.5*(2*l+1)*fak(l-m)/fak(l+m);
end

function P = legendre(l,m)
	P = 0;
	error
end

function L = laguerre(a, n, x)
	L = (x^(-a)*e^x)/fak(n) * dn_dxn * (e^(-x)*x^(n+a));
	error
end

function F = fak(N)
	if (N > 1)
		F = N*fak(N-1);
	else
		F = 1.0;
	end
end

