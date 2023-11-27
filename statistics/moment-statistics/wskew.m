% Do 11. Feb 19:33:24 CET 2016
% Karl Kastner, Berlin
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
%
%% skewness of a weighted set of samples
%% function sk = wskew(w,x)
function sk = wskew(w,x)
	w  = w./sum(w);
	mu = wmean(w,x);
	sd = wstd(w,x);
	c3 = w'*(x-mu).^3;
	sk = c3./(sd.^3);
end

