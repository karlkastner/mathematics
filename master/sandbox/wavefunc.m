% 2011 Oct  7 13:28 MSK
% Karl Kästner, Berlin

function phi = wf(r)
	a0 = 1;
	%phi = 1/sqrt(pi)*a0^(-3/2)*exp(-abs(r)/a0);
	phi = (4*pi*r.^2).*(1/sqrt(pi)*a0^(-3/2)*exp(-abs(r)/a0)).^2;
end % wf

