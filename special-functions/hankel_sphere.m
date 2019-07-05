% Fri Jun 28 08:09:37 UTC 2013
% Karl Kästner, Berlin
%
%% spherical Hankel function for the far field (incident plane wave)
%% first kind
function h = hankel_sphere(kR,m)
	h = bessel_sphere(kR, m) + 1i*neumann_sphere(kR,m);
end % function hankel_sphere

