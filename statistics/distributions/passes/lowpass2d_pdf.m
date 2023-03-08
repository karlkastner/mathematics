% Fri 22 Apr 14:05:05 CEST 2022
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
% unnormalized (radial) density of the pth-order lowpass in two dimensions
% continuous space
function S = lowpass2d_pdf(fr,a,order)
	S = zeros(size(fr));
	for idx=1:numel(fr)
		S(idx) = integral(@(r) besselj(0,2*pi*abs(r*fr(idx))).*r.*exp(-a*abs(r)),0,inf);
	end
	% S0 = int r S J0(r k) dk = 1/a^2
	S0 = 1./a^2;
	% scale S(0) to 1
	S = S/S0;
	if (nargin()>2 && ~isempty(order))
		S = S.^order;
	end
	% TODO normalize
end

