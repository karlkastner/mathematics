% Do 3. Mär 21:06:44 CET 2016
% Karl Kastner, Berlin
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
%%
%% probability density of and unknown distribution
%% with mean mu, standard deviation sigma, and third and fourth cumulants
%% c.f. Rao 2010
function f = edgeworth_pdf(x,mu,s,sk,ku,n)
	% normalise x
	x = (x-mu)/s;
	H = @(n,x) hermiteH(n,x);
	f = normpdf(x).*(  1 ... 
			 + 1/6*sk*H(3,x) ...
			 + 1/24*(ku-3)*H(4,x) ...
			 + 1/72*sk^2*H(6,x));
end

