% Thu  6 Jan 15:53:31 CET 2022
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
%% normaliztation scale of the spatial bandpass density
%
% function [Sc] = bandpass1d_continuous_pdf_scale(fc,p,pp,numeric)
function [Sc] = bandpass1d_continuous_pdf_scale(fc,p,pp,numeric)
	if (nargin()<3)
		pp = [];
	end
	if (nargin()<4)
		numeric = false;
	end
%	else
	if (isempty(pp) && ~numeric)
		kc = 2*pi*fc;
		% IS = kc./(2*pi) * 4.^(1-p)./(2*p-1).*pi.*p.*binom(2*p-1,p-1)
		IS = kc.*2.^(1-2*p).*binomial(2*(p-1),p-1);
		if (isnan(IS))
			% stirlings approximation
			IS = kc./sqrt(4*pi*(p-1));
		end
		if (0 == IS)
			Sc = realmax;
		else
			Sc = 1./IS;
		end
	else
		% avoid overflow
		% TODO these limits are for q = 1
		tol = 1e-5^(1/(4*p));
		fl = fc*(1-sqrt(1 - tol^2))/tol;
		fr = fc*(1+sqrt(1 - tol^2))/tol;
		Sc = 1./quad(@(fx) bandpass1d_continuous_pdf(fx,fc,p,0,pp),fl,fr);
	end
	if (~issym(p))
		Sc(p<0.5) = NaN;
	end
end

