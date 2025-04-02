% Fri 22 Apr 14:05:05 CEST 2022
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
%
%% function S = lowpass2dpdf(fx,fy,fc,order,normalize)
%%
%% spectral density of the two-dimensional lowpass filter
%% in continuous space
%%
%% ouput :
%% S : spectral density
%%
%% input :
%% fx, fy : frequencies at which to evaluate the filter
%% fc     : characteristic frequency
%% p      : order
%%
%% note that for p = 3/4, the autocorrelation function becomes exponential
%%
%% R = exp(-a*sqrt(x^2 + y^2))
%% c.f. Baddour and table of Wolfram Hankel transform:
%%
%% S(f) = int R(r) J0(r 2 pi f) r dr 
%%      = int exp(-a r) J0(2 pi f r)
%%      = 2 pi a^2 (a^2 + 4*pi^2*f^2)^(-3/2)
function S = lowpass2dpdf(fx,fy,fc,p,normalize)
	if (nargin()<4)
		p = 1;
	end
	if (nargin()<5)
		normalize = 0;
	end
	fr = hypot(fx,fy);
	S  = lowpass1dpdf(fr,fc,p,false);
	switch (normalize)
	case {0}
		% do not normalize
	case {1}
		% before normalization S(fc) = 1,
		% after normalization S(fc) = Sc,
		% so the normalization factor is Sc = max(S)
		Sc = lowpass2dpdf_max(fc,p);
		S = Sc.*S;
	end
end % lowpass2dpdf

