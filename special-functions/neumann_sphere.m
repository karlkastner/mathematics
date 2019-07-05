% Fri Jun 28 08:10:05 UTC 2013
% Karl Kästner, Berlin
%
%% spherical Neumann function
%% Bessel function of the second kind
function n = neumann_sphere(kR,m)
	% far field (incident plane wave)
	%n = 1/kR*sin(kR - 0.5*(m+1)*pi);
	Y = bessely(m+0.5,kR);
	n = sqrt(0.5*pi/kR)*Y;
end


