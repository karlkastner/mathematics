% Mon 19 Dec 17:03:02 CET 2022
% Karl Kästner, Berlin
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
% function [w, dof] = gausswin2(n, nf)
%
% Gaussian filter window in 2D with sum of weights normalized to 1
% input :
% 	n  = [nx,ny]   : size of output window
% 	nf = [nfx,nfy] : cut off radius of Gaussian
% 			 the cut off radius can be non-integer
% output:
%	w  (nx * ny)   : window matrix
%	dof            : estimated degrees of freedom,
%			 in case samples smoothed by the window are independent
%
% TODO allow for rotation in case of anisotropy
function [w, dof] = gausswin2(n, nf)
	[fx, fy, fr, ft] = fourier_axis_2d([1,1],n);

	% Gaussian window
	if (1 == length(nf))
		w = normpdf(fr,0,nf/log(4));
	else
		wx = normpdf(fr,nf(1)/log(4));
		wy = normpdf(fy,nf(2)/log(4));
		w  = cvec(wx)*rvec(wy);
	end

	% normalize sum of weigths to 1
	w = w/sum(w(:));

	% degrees of freedom of window
	% when values are uncorrelated
	if (nargout()>1)
		dof = sum(w(:)).^2./(sum(w(:).^2));
	end
end

