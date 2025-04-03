% Tue  4 Feb 16:33:31 CET 2025
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
function sd = binormpdf_std(p,mu1,mu2,sd1,sd2)
	% \sigma^2 = p sigma1^2 + (1-p) sigma2^2 + p*(1-p)*(mu0 - mu1)^2
	s2 = binormpdf_var(p,mu1,mu2,sd1,sd2);
	sd = sqrt(s2);
end

