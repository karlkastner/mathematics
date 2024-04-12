% Fri 22 Apr 13:28:53 CEST 2022
% Karl Kästner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% not normalized, max (S) = 1;
% function Sb = bandpass2d_continous_pdf_hankel(L,n,a,order,varargin)
function Sb = bandpass2d_continous_pdf_hankel(L,n,a,order,varargin)
	% lowpass density
	Sl = lowpass2d_continuous_pdf_hankel(L,n,a,[],varargin{:});
	% bandpass density
	Sb = 4*Sl.*(1.0-Sl);
	% higher order
	if (nargin()>2 && ~isempty(order))
		Sb = Sb.^order;
	end
	% normalize
	df = 1/L;
	Sb = 2*Sb./(sum(Sb)*df);
end

