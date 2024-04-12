% Wed 22 Nov 09:36:22 CET 2023
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function D = derivative_matrix_k_1d(n,L,k)
	dx = L/(n-1);
	s = 1/dx^k;
	l = 2*fix((k+1)/2)+1;
	d = s*derivative_matrix_kernel(l,k,0);
	D = spdiags(ones(n,1)*d,-(l-1)/2:(l-1)/2,n,n);
	
	% extrapolation at boundary, constant for even derivatives, linear for odd derivatives
	% linear extrapolation requires longer kernel
	for idx=1:(l-1)/2
		d = s*derivative_matrix_kernel(l,k,-(l+1)/2+idx);
		D(idx,1:l) = d;
		d = s*derivative_matrix_kernel(l,k, (l+1)/2-idx);
		D(end-idx+1,(end-l+1):end) = d;
	end
end

