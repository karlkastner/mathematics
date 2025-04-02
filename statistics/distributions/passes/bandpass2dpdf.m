% Fri 22 Apr 13:28:53 CEST 2022
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
%% function Sb = bandpass2dpdf(fr,a,order,normalize)
%%
%% ouput:
%% S : spectral density of the bandpass in two-dimensions
%% 
%% input:
%% fx, fy : frequencies where to evaluate the spectral density
%% fc : characteristic frequency
%% p  : order
%% normalize : if 0 : S(fc) = 1
%              if 1 : int int S dfx dfy = 1
function Sb = bandpass2dpdf(fx,fy,fc,order,normalize)
	% lowpass density
	Sl = lowpass2dpdf(fx,fy,fc,1,0);
	% first order bandpass density
	Sb = 4*Sl.*(1.0-Sl);
	% higher order
	if (nargin()>2 && ~isempty(order))
		Sb = Sb.^order;
	else
		order = 1;
	end
	% normalization
	switch (normalize)
	case {1}
		% without normalization S(fc) = 1
		% with normalization S(fc) = Sc = max(S)
		% so the normalization factor is Sc
		Sc = bandpass2dpdf_max(fc,order);
		Sb = Sc.*Sb;
	end
end

