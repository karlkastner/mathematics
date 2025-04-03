% 2025-03-17 20:35:42.113447194
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
function [a,b,p] = generalized_gamma_cmoment2par(mu,s2,sk)
	% initial condition
	[a,b] = gamma_moment2par(mu,s2);
	par0 = [a,b,1];
	par = lsqnonlin(@(par) [generalized_gamma_mean(p(1),p(2),p(3)) - mu,
			  generalized_gamma_var(p(1),p(2),p(3)) - s2,
			  generalized_gamma_skewnewss(p(1),p(2),p(3)) - sk], ...
			  par0 ...
		 );
	a = par(1);
	b = par(2);
	p = par(3);
end

