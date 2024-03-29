% Mon 19 Dec 17:03:02 CET 2022
% Karl Kästner, Berlin
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
% gaussian window in 1D
function [w, dof] = gausswin1(n, nf)
	fx = fourier_axis(1,n);

	% Gaussian window
	% 0.5 = exp(-1/2*(nf/2)^2/s^2)
	% -2*log(0.5) = (nf/2)^2/s^2
	% -8*log(0.5) = nf^2/s^2
	% s^2 = -nf^2/(8*log(0.5))
	% s^2 =  nf^2/(8*log(2))
	% s   =  nf/sqrt(8*log(2))
	s = nf/sqrt(8*log(2));
	w = normpdf(fx,0,s);

	% normalize sum of weigths to 1
	w = w/sum(w(:));

	% degrees of freedom of window
	% when values are uncorrelated
	dof = sum(w(:)).^2./(sum(w(:).^2));
end

