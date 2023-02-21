% Sat 28 Jan 14:58:19 CET 2023
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
%
% n.b. this has an error of about 12% for nmsk -> n, as the discrete transform
% of a rectangular pulse slightly deviates from the sinc-function
%
% function [Leff,S_bmsk] = effective_mask_size(bmsk,L,angle_deg)
function [Leff,S_bmsk] = effective_mask_size(bmsk,L,angle_deg)
	n       = size(bmsk);
	% sinc(scale)^2 = 0.5
	scale   = 0.443;
	[fx,fy] = fourier_axis_2d(L,n);
	S_bmsk  = abs(fft2(bmsk)).^2;
	[Sr,fr] = periodogram_radial(S_bmsk,L);
	Sr      = Sr.normalized;
	S_bmsk  = fft_rotate(S_bmsk,angle_deg);

	% find place, where the density drops to 0.5
	fdx    = monotoneous_indices(S_bmsk(:,1),'descending');
	fx05   = interp1(S_bmsk(fdx,1),fx(fdx),0.5*S_bmsk(1,1),'linear');
	fx05   = fx05;
	Leff.x = scale./fx05;

	fdx    = monotoneous_indices(S_bmsk(1,:),'descending');
	fy05   = interp1(S_bmsk(1,fdx),fy(fdx),0.5*S_bmsk(1,1),'linear');
	fy05   = fy05;
	Leff.y = scale./fy05;

	fdx    = monotoneous_indices(Sr,'descending');
	fr05   = interp1(Sr(fdx),fr(fdx),0.5*Sr(1),'linear');
	fr05   = fr05;
	Leff.r = scale./fr05;
end

