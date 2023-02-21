% Tue 21 Feb 14:16:29 CET 2023
% Karl Kastner, Berlin
%
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
% demonstration, that averaging densities with the same distribution
% but different regularity results in a density that is more pointed
% and has heavier tales than the underlying distribution
%
%
%% generalized binomial coefficient, working for non-integer arguments,
%% in contrast to the matlab buildin function nchoosek
function b = binomial(n,k)
	if (n+1 < 0 || k+1<0 | n-k+1<0)
		b = NaN;
	else
		b = exp(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1));
	end
end

