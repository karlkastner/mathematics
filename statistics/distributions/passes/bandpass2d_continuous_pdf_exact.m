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
%% function Sb = bandpass2d_continuous_pdf(fr,a,order)
%% not normalized, max (S) = 1;
function Sb = bandpass2d_continuous_pdf(fr,a,order)
	% lowpass density
	Sl = lowpass2d_continuous_pdf(fr,a,[]);
	% bandpass density
	Sb = 4*Sl.*(1.0-Sl);
	% round of error
	Sb = max(0,Sb);
	% higher order
	if (nargin()>2 && ~isempty(order))
		Sb = Sb.^order;
	end
end

