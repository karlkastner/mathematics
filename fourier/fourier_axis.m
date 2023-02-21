% Di 12. Jan 09:26:57 CET 2016
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
%% return axis of frequencies and periods for the discrete fourier transform
%% as computed by fft (matlab-style)
%%
%% input:
%% X : sample locations (equal interval)
%% L : length of samples
%% n : number of samples
%%
%% output :
%% f    : frequencies
%% T    : periods
%% mask : mask for 1/2 of the fourier transform
%%        (as both halves are complex conjugates)
%% N    : frequency id
% function [f T mask N] = fourier_axis(X)
% function [f T mask N] = fourier_axis(L,n)
%
% note that fa(fa(x)) = x
function [f, T, mask, N] = fourier_axis(varargin)
	if (1 == nargin())
		X = varargin{1};
		n = length(X);
		L = (max(X)-min(X))*n/(n-1);
		%L = (max(X)-min(X))*n/(n-1);
	else
		L = varargin{1};
		n = varargin{2};
	end
	T = zeros(n,1);
	k = ceil(n/2);
	mask = false(n,1);
	% even
	if (0 == mod(n,2))
		%N = [0:n/2,-(n/2-1:-1:1)]';
		%mask(1:n/2+1) = true;
		N = [0:n/2-1,-(n/2:-1:1)]';
		mask(1:n/2) = true;
	else
		N = [0:(n-1)/2,-((n-1)/2:-1:1)]';
		mask(1:(n+1)/2) = true;
	end

	% frequency
	f = (1/L)*N;
	% wave length
	T = 1./f;

	% variance of the fourier coefficients: (coefficients normalised by n)
	% c.f: fourier-analysis-for-beginners-ch8-statistical-description-of-fourier-coefficients-thibos-2003
	% s2 = sigma(Y);
	% s2f(1)     = 2*s2/n;
	% s2f(2:end) = 2*s2/n;
	
end

