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
%% radial spectral density of the pth-order lowpass in two dimensions
%% and continuous space
%% with autocorrelation
%% R = exp(-a*sqrt(x^2 + y^2))
%%
%% determined by numerical integration of the exact expression
%%
% function S = lowpass2d_pdf(fr,a,order)
%
function Sr = lowpass2d_continuous_pdf(fr,a,order,normalize)
	Sr = zeros(size(fr));
	for idx=1:numel(fr)
		Sr(idx) = integral(@(r) besselj(0,(2*pi)*r*fr(idx)).*r.*exp(-a*r),0,inf);
	end
	% S0 = int r S J0(r k) dk = 1/a^2
	S0 = a.^2;
	% scale S(0) to 1
	Sr = Sr/S0;
	Sr=Sr/Sr(1);
	if (nargin()>2 && ~isempty(order))
		Sr = Sr.^order;
	end
	if(nargin()>3 && normalize)
		% TODO S(0) should be weighted by 0.5
		Sr = Sr/(sum(Sr)*(fr(2)-fr(1)));
	end
end

