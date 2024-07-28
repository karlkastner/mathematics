% 2024-06-30 11:55:28.738845211 +0200
function [fc,Sc] = cauchywrappedpdf_mode(f0,s)
	r = 2*abs(f0)*sqrt(f0.^2 + s.^2) - f0.^2 - s.^2;
	if (r < 0)
% (- 2*f0*(f0^2 + s^2)^(1/2) - f0^2 - s^2)^(1/2)
%   (2*f0*(f0^2 + s^2)^(1/2) - f0^2 - s^2)^(1/2)
%-(- 2*f0*(f0^2 + s^2)^(1/2) - f0^2 - s^2)^(1/2)
%  -(2*f0*(f0^2 + s^2)^(1/2) - f0^2 - s^2)^(1/2)
%		2*abs(f0)*sqrt(f0.^2 + s.^2) 
%		- f0.^2 
%		- s.^2
%		pause
		fc = 0;
	else
		fc = sqrt(r);
	end
	Sc = cauchywrappedpdf(fc,f0,s);
end

