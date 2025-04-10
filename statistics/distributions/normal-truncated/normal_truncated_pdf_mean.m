% Tue 11 Feb 15:32:51 CET 2025
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
function mean_ = normal_truncated_pdf_mean(mu,sd,lb)
	error('broken')
	% shift lb to 0
	mu = mu-lb;
	h = -mu/sd;
	r = normpdf(h)/(1-normcdf(h));
	mean_ = mu + r*sd;
end

