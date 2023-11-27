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
% spectral density of the discrete phase noise integration process
% (damped oscillator with phase noise)
% 	(z + a1 dz/dx + a2 d^2z/dx^2) == e
function [S1,AA] = damped_oscillator_2d_discrete_pdf(n,L,a1,a2,a3)
	Ix  = speye(n(1));
	Iy  = speye(n(2));
	D1x = derivative_matrix_1_1d(n(1),L(1),2,'circular','circular');	
	D2x = derivative_matrix_2_1d(n(1),L(1),2,'circular','circular');
	D2y = derivative_matrix_2_1d(n(2),L(2),2,'circular','circular');
	D1x = kron(D1x,Iy);
	D2x = kron(D2x,Iy);
	D2y = kron(Ix,D2y);
	% TODO compute D4 directly
	D4x = D2x*D2x;
	I   = speye(prod(n));
%	AA = A'A
%          = (I + a1*D1x + a2*D2x + a3*D2y)*(I - a1*D1x + a2*D2x + a3*D2y);
%          = (       I -    a1*D1x +    a2*D2x + a3*D2y 
%             + a1 D1x -  a1^2 D1x + a1 a2 D3x + a1 D1xD2y
%             + a2 D2x - a1 a2 D3x +  a2^2 D4x + a2 a3 D2xD2y
%             + a3 D2y - a1 a3 D1x D2y + a2 a3 D2xD2y + a3^2*D4y)
%       A = (I + (2*a2 - a1^2)*D2x + a2^2*D4x - a3*D2y);
	A = (I + a1*D1x + a2*D2x + a3*D2y);
	AA = A'*A;
	z = zeros(prod(n),1);
	z(1,1) = 1;
	if (1)
		acf2 = bicgstabl(AA,z);
		%acf2 = minres(AA,z);
	else
		acf2 = AA \ z;
	end
	acf2 = reshape(z,n);
	S2   = ifft2(acf2);
	S1   = sqrt(S2);
	df   = 1./L;
	S1   = 2*S1/(sum(S1,'all')*df(1)*df(2));
end
