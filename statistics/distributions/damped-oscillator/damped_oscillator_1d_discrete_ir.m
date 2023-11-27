% Wed 25 Oct 11:04:13 CEST 2023
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
function [ir,A] = damped_oscillator_1d_discrete_ir(n,L,a1,a2)
	I  = speye(n);
	D1 = derivative_matrix_1_1d(n,L,2,'circular','circular');	
	D2 = derivative_matrix_2_1d(n,L,2,'circular','circular');
	A = (I + a1*D1 + a2*D2);
%	A = (I + a1*D1 + a2*D2)*(I - a1*D1 + a2*D2);
%	A = (I - a1*D1 + a2*D2 + a1 D1 - a1^2 D1^2 + a1 a2 D1 D2 + a2 D2 - a1 a2 D1 D2 + a2^2 D2^2)
%	D4 = D2*D2;
%	D2_ = D1*D1;
%	A_ = (I + (2*a2*D2 - a1^2*D2_) + a2^2*D4)
%	A_ = (I + (2*a2 - a1^2)*D2 + a2^2*D4)
%	sum(sum(abs(A-A_)))
%	A = (I + a1*D1)*(I - a1*D1)*(I + a2*D2);
	z = zeros(n,1);
	z(1) = L;
	% impulse response
	ir = A \ z;
%	ir=acf/acf(1);
end
