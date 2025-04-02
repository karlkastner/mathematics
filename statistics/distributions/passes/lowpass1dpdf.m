% Sat 26 Jun 21:04:19 CEST 2021
% Karl Kastner, Berlin
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
% function [S_lp, S_lp1] = lowpass1d_continuous_pdf(fx,fc,p,normalize)
%
%% function [S_lp, S_lp1] = lowpass1dpdf(fx,f0,p,normalize)
%%
%% spectral density of the one-dimensional spatial, i.e. two-sided, low-pass filter
%% identical to the radial density of the two-dimensional low-pass-filter
%%
%% fx : frequencies at which to evaluate density
%%
%% f0 : frequency scale
%% p  : order of the filter
%% normalize : if 1, area under the density between 0 and infty is 1,
%%	      otherwise, S(0) = 1
function [S_lp, S_lp1] = lowpass1dpdf(fx,f0,p,normalize)
	if (nargin()<3)
		p = 1;
	end
	if (nargin()<4)
		normalize = true;
	end
	S_lp1 = 1./(1 + fx.^2./f0.^2).^2;
	S_lp  = S_lp1.^(p);
	switch (normalize)
	case {0}
		% no normalization
	otherwise
		% note that before nomralization S(0) = 1 and after normalization
		% S(0) = max(S), so the normalization factor is max(S)
		Sc = lowpass1dpdf_max(f0, p);
		S_lp = Sc.*S_lp;
	end
end

