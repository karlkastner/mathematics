% 2025-02-27 19:13:01.558542485 +0100
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
% 
%% determine correlation length of the process
%%	z = exp(lz)
%%	z  ~ lognormal(mu,sd)
%%	lz ~ normal(lmu,lsd)
%%	R_lz(x,y) = exp(-sqrt(x^2 + y^2)/lz)
%% correlation lengths
%%	Rz(t)   = exp(-1)
%%	Rzl(tl) = exp(-1)
%	
% TODO use s not ls
function lt = geometric_ou_correlation_length_lin2log(ls,t)
	lt = -t/log((log(exp(ls^2) + exp(1) - 1) - 1)/ls^2);
end

