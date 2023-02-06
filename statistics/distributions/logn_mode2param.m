% 2022-03-30 12:46:31.571547062 +0200

function [lmu,lsd] = logn_mode2param(xm,ym)
	% sqrt(2*pi)*lc*Sc = exp(-1/2*s^2)/s
	lsd2 = lambertw(1./(2*pi*xm.^2.*ym.^2));
	lsd = sqrt(lsd2);
	lmu = log(xm)+lsd2;
end

