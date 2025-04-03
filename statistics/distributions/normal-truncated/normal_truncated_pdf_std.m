% Tue 11 Feb 15:33:00 CET 2025
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
function sd = normal_truncated_pdf_std(mu,sd,lb,ub)
	mu = mu-lb;
	h = -mu/sd;

	r = normpdf(h)/(1-normcdf(h));
	% Z = normcdf(a)
	% s2  = 1 - a*normpdf(a)/Z - (normpdf(a)/Z).^2
	s2 = sd^2*(1 + r*(h - r)); 
	sd = sqrt(s2);
end

