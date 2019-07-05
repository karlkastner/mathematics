% Fri Jun 28 08:10:29 UTC 2013
% Karl Kästner, Berlin
%
%% spherical Bessel function of the first kind
function j = bessel_sphere(kR,m)
	% far field (incident plane wave)
	% j = 1/kR*cos(kR - 0.5*(m+1)*pi);
	J = besselj(m+0.5,kR);
	j = sqrt(0.5*pi/kR)*J;
end % function j

