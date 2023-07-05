% Sat 26 Jun 21:04:19 CEST 2021
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
% function [S_lp, S_lp1] = lowpass1d_continuous_pdf(fx,fc,p,normalize)
function [S_lp, S_lp1] = lowpass1d_continuous_pdf(fx,fc,p,normalize)
	if (nargin()<3)
		p = 1;
	end
	if (nargin()<4)
		normalize = true;
	end
	S_lp1 = 1./(1 + fx.^2./fc.^2);
	S_lp  = S_lp1.^p;
	switch (normalize)
	case {0}
		% no normalization
		Sc = 1;
	otherwise
		Sc = lowpass1d_continuous_pdf_scale(fc, p);
	end
	S_lp = Sc.*S_lp;
end

