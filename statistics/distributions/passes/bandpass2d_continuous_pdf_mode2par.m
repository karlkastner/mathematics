% Fri  7 Jan 16:11:57 CET 2022
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
%
%% transform mode (maxima) of the bandpass spectral density into the paramter
%% of the underlying distribution 
%
% function [p] = spectral_density_bandpass_max2par(fc,Sc,p0)
function [p1,p2] = bandpass2d_continuous_pdf_mode2par(fc,Sc,par0,L,n) %par0,pp)
	if (nargin()<3||isempty(par0))
		par0 = [6*fc,1];
	end
	if (nargin()<4)
		L = 100/fc;
	end
	if (nargin()<5)
		n = round((L*fc)^2);
	end
%	if (nargin()<4)
%		pp = [];
%	end

	%n = 20;
	%lc = 1/fc;
	%L = n*lc;
	%df = 1/L;
	%fx = fourier_axis(L,n)';

	opt = struct();
%	opt.Algorithm = 'levenberg-marquardt';
	par = lsqnonlin(@(par) fun(par) - [fc,Sc], par0,[0,0],[],opt);
%	par = lsqnonlin(@(par) fun(par) - [fc,Sc], par,[0,0],[],opt);
	p1 = par(1);
	p2 = par(2);

	function mode = fun(par)
		[mode(1),mode(2)] = bandpass2d_continuous_pdf_mode(par(1),par(2),L,n);
	end
end

