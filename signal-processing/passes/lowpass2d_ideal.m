% Wed 27 Apr 10:55:29 CEST 2022
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
%% lowpass filter the input x in the Frequency Domain
%%
%% TODO no need to provide dx, follows from size of x
%% function [y,S,R,r]=lowpass2d_ideal(x,L,dx,varargin)
function [y,S,R,r]=lowpass2d_ideal(x,L,a,varargin)
	n = size(x);
	[fx,fy,fr] = fourier_axis_2d(L,n);
	%[S,R,r] = lowpass2d_pdf(size(x),dx,L,varargin{:});
	[S] = lowpass1d_continuous_pdf(fr,a,varargin{:});
	y = ifft2(sqrt(S).*fft2(x));
end

