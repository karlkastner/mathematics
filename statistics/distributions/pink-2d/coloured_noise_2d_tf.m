% 2024-01-22 10:07:14.078606190 +0100
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%% transfer function of coloured noise
function [T,fx,fy,frr] = coloured_noise_2d_tf(n,L,p);
	% spectral density of coloured noise
	[fx,fy,frr] = fourier_axis_2d(n,L);
	% S = 1./frr^p;
	% T = sqrt(S);
	% transfer function
	T = frr.^(-0.5*p);
	T(1,1) = 0;
end

