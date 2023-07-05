% Do 11. Feb 18:11:02 CET 2016
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
%% edgeworth expansion of an unknown cumulative distribution
%% with mean mu, standard deviation sigma, and third and fourth cumulants
%% c.f. Rao 2010
function F = edgeworth_cdf(mu,sigma,c3,c4,x)
	% standardise
	x  = (x-mu)/sigma;
	p1 = -1/6*c3*(x.^2-1);
	% these are hermite polynomials
	p2 = -x.*(1/24*c4*(x.^2-3) ...
	     + 1/72*c3^2*(x.^4 - 10*x.^2 + 15));
	F  = normcdf(x) + normpdf(x).*(p1 + p2);
end

