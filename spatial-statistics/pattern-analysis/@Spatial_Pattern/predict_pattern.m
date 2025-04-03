% 2025-03-26 10:56:39.606912627 +0100
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
function predict_pattern(obj)
	T = interp1(obj.f.r,obj.T.radial,obj.f.rr,'linear',0);
	obj.b_.lin = ifft2(T.*fft2(obj.source));
	% strip spurious imaginary part introduced by finite machine precision
	obj.b_.lin = real(obj.b_.lin);
	% threshhold
	q = quantile(obj.b_.lin,1-obj.stat.coverage,'all');
	obj.b_.lin_thresh = obj.b_.lin > q;
	% goodness of fit
	obj.stat.fit.b_.lin.r2        = corr(obj.b_.lin(:),obj.b_.square(:)).^2;
	obj.stat.fit.b_.lin_thresh.r2 = sign_to_pearson(corr(obj.b_.lin_thresh(:),obj.b_.thresh(:))).^2;
end

