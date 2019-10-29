% Thu Jan 19 03:21:47 MSK 2012
% Karl Kästner
% 
% wave function of the unconfined hydrogen atom

function [Psi RR TT PP X Y Z] = hydrogen_wf(n,l,m, r, theta, phi)
	x = linspace(-25,25,100);
	y = linspace(-25,25,100);
	z = linspace(-50,0,100);
	[X Y Z] = meshgrid(x,y,z);
	r = sqrt(X.^2 + Y.^2 + Z.^2);
	phi   = atan2(Y,X);
	theta = acos(Z./r); 
	[theta phi r] = cart2sph(X,Y,Z);
	theta=pi/2-theta;
	RR = r;
	PP = phi;
	TT = theta;

	  psi_r = radial(n, l, r);
	  % y = psi_c * psi_a = spherical harmonics
	  psi_c = colatitude(l, m, theta);
	  psi_a = azimuthal(m, phi);
	  Psi = psi_r.*psi_c.*psi_a;
	  %Psi = sqrt(conj(Psi).*Psi);
	  Psi = (conj(Psi).*Psi);
	return;

	  Psi = kron(kron(psi_r, psi_c), psi_a);
	  RR  = kron(kron(r,ones(size(theta))), ones(size(phi)));
	  TT  = kron(kron(ones(size(r)),theta), ones(size(phi)));
	  PP  = kron(kron(ones(size(r)),ones(size(theta))), phi);

	  %X = kron(kron(r, sin(theta)), cos(phi));
	  %Y = kron(kron(r, sin(theta)), sin(phi));
	  %Z = kron(kron(r, cos(theta)), ones(size(phi)));
	  X = kron(kron(r, cos(phi)), sin(theta));
	  Y = kron(kron(r, sin(phi)), sin(theta));
	  Z = kron(kron(r, ones(size(phi))), cos(theta));

	  Psi = reshape(Psi,length(r),length(theta),length(phi));
	  X  = reshape(X,length(r),length(theta),length(phi));
	  Y  = reshape(Y,length(r),length(theta),length(phi));
	  Z  = reshape(Z,length(r),length(theta),length(phi));
end % hydrogen wf

function f = radial(n, l, r)
	a0 = 1;
	rho = 2*r/(n*a0);
	% cubic zu viel on wiki formula
	f = sqrt( (2/(n*a0))^3 * factorial(n - l - 1)/(2*n*factorial(n+l)) ) ...
	    .* exp(-rho/2).*rho.^l.*assoc_laguerre(n-l-1,2*l+1,rho);
	    %.* exp(-rho/2).*rho.^l.*assoc_laguerre(2*l+1,n-l-1,rho);
end % function radial

function f = colatitude(l, m, theta)
	f =   sqrt((2*l+1)/2 * factorial(l-m)/factorial(l+m)) ...
	     * assoc_legendre(l,m,cos(theta)) ...
	    .* (1-cos(theta).^2).^(abs(m)/2); % unsure about this line
end % function colatitude

function f = azimuthal(m, phi)
	f = 1/sqrt(2*pi) * exp((1i*m).*phi);
end % function azimuthal

