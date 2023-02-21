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
% function [w, nf2] = circwin(n, nf)
%
% circular (elliptical) window in 2D
%
% input:
%	n = [nx,ny] : size of image containing the window
%	nf = nr or [nfx,nfy] : windows size
%		radius of circle or axis of ellipse that is non-zero
%		allowing for non-integer radii
%
% output :
%	w   : window with sum of weights normalized to 1
%	nf2 : number of non-zero bins in window
%
% TODO allow for rotation in case of anisotropy
function [w, dof] = circwin(n, nf)
	[fx, fy, fr, ft] = fourier_axis_2d([1,1],n);

	if (1 == length(nf))
		nf = [nf,nf];
	end

	% integer part
	nfi = floor(nf);
	p   = nf-nfi;

	% circular window
	% note : the area seems to be only first order accurate because the
	%	 disk of the circle is discontinuous
	[wi,pi] = int_1d_gauss_2();
	%[wi,pi] = int_1d_nc_2();
	di = pi*[-0.5;0.5];
	w = zeros(n);
	fx = cvec(fx);
	fy = rvec(fy);
	for idx=1:length(di)
	for jdx=1:length(di)
		w = w + wi(idx)*wi(jdx)*(hypot((fx+di(idx))/nf(1),(fy+di(jdx))/nf(2)) <= 1);
	end
	end

	% normalize sum of weights to 1
	w = w./sum(w,'all');

	% degrees of freedon, when samples are independent
	dof = 1./sum(w.*w,'all');
end % circwin

