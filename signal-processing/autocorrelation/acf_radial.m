% 2022-07-11 09:54:04.582898046 +0200
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
% function R = radial_acf(R,L)
function R = radial_acf(R,L)
	n = size(R);
	if (nargin()<2)
		L = n;
	end
	fx = fourier_axis(L(1),n(1));
	fy = fourier_axis(L(2),n(2));
	fmax = max(max(fx),max(fy)); 
	df   = 1./L;
	dfi  = sqrt(df(1).*df(2));
	fr   = hypot(fx,fy');
	frmax = hypot(max(fx),max(fy));
	fri  = (0:dfi:ceil(frmax))';

	% this is a linear interpolation
	% each value is proportionally split between the next lower and next larger bin
	% integer part
	rat = fr/dfi;
	bin = floor(rat);
	p   = 1-(rat-bin);
	% 0 is bin 1
	bin = bin+1;

	s = [length(fri)+2,1];
	sumR = ( accumarray(bin(:),p(:).*R(:),s,@sum) ...
	       + accumarray(bin(:)+1,(1-p(:)).*R(:), s, @sum) );
	% number of bins
	count = ( accumarray(bin(:),p(:),s,@sum) ...
	       + accumarray(bin(:)+1,(1-p(:)), s, @sum) );

	R = sumR./count;
	R = R/R(1);
end
