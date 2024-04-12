% Wed  1 Dec 18:55:18 CET 2021
% Karl Kästner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% spectral density of a fourier series where the phase undergoes brownian motion
%% with standard deviation s per unit distance
%
%function [S,I] = phase_drift_pdf(fx,f0,s,normalize)
function [S,I] = phase_drift_pdf(fx,f0,s,normalize)
	if (nargin()<4)
		normalize = 0;
	end
	if (issym(s))
		pi_ = sym(pi);
	else
		pi_ = pi;
	end
	p = pi_*s.^2;
	% this is already analytically normalized
	S = 4*p.*(p.^2 + fx.^2./f0.^2 + 1)./(2*pi_*f0.*(4*p.^2 + (p.^2 + fx.^2./f0.^2 - 1).^2));


	switch (normalize)
	case {1} % nuermically
		% TODO, use quad here
		I = spectral_density_area(fx,S);
		S = S./I;
	end
end % spectral_density_brownian_phase

