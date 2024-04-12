% Fri  7 Jan 12:44:47 CET 2022
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
%% transform mode to parameters of the brownian phase spectral density
%
% function [f0, s] = phase_drift_pdf_mode2par(fc,Sc)
function [f0, s] = phase_drift_pdf_mode2par(fc,Sc,numerical)
	if (nargin()<3 || ~numerical)
		%r = root(z^3 - z^2 + z*(4*Sc^2*fc^2*pi^2 - 2) - 4*Sc^2*fc^2*pi^2, z, 1)
		Sfp2 = Sc*Sc*fc*fc*pi*pi;
		rp = [1, -1, 4*Sfp2 - 2, -4*Sfp2];
		r = roots3(rp);
		% choose real root
		r	
		r = r(1);
		% by round off error, r can sometimes ha
		r = real(r);
		den1 =  2*( (2*Sfp2 - 1)*r.*r ...
	                    + r ...
		            + Sfp2*(-6 + 8*Sfp2) ...
		          );
		s = 2*sqrt(Sc.*fc./sqrt(den1));
		den2 = (2*sqrt(1 + pi^2*s.^4) - s.^4*pi^2 - 1);
		f0 = fc./sqrt(den2);
		if (den1  < 0 || den2 < 0)
			[f0, s] = phase_drift_pdf_mode2par(fc,Sc,true);
		end
	else
		fc = double(fc);
		Sc = double(Sc);
		par = ([fc,0.5]);
		par = lsqnonlin(@objfun,[par(1),par(2)],[0,0]);
		f0 = par(1);
		s = par(2);
	end
function res = objfun(par)
	[fc_,Sc_] = phase_drift_pdf_mode(par(1),par(2));
	res = [fc_-fc,Sc_-Sc];
end
end % phase_drift_pdf_mode2par

