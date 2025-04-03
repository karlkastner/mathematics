% 2025-02-19 16:55:56.692660851 +0100
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
function r = corr_angle(a,b)
	sa = sin(a);
	ca = cos(a);
	mua = atan2(mean(sa),mean(ca));
	sb = sin(b);
	cb = cos(b);
	mub = atan2(mean(sb),mean(cb));

	sa_ = sin(a-mua);
	sb_ = sin(b-mub);
	
	r = mean(sa_.*sb_)./sqrt(mean(sa_.*sa_.*sb_.*sb_));
end
