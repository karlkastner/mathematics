% Thu 26 Oct 15:14:11 CEST 2023
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
% spectral density of the phase noise integration process
% note that the density and acf are real and symmetric,
% while the ir is asymmetric and the tf complex
function [S1,AA] = damped_oscillator_1d_discrete_pdf(n,L,a1,a2)
	I  = speye(n);
	D1 = derivative_matrix_1_1d(n,L,2,'circular','circular');	
	D2 = derivative_matrix_2_1d(n,L,2,'circular','circular');
	D4 = D2*D2;
	% note :
	% AA = A'*A
        %    = (I + a1*D1 + a2*D2)*(I - a1*D1 + a2*D2)
        %    = (I - a1*D1 + a2*D2 + a1 D1 - a1^2 D1^2 + a1 a2 D1 D2 + a2 D2 - a1 a2 D1 D2 + a2^2 D2^2)
	AA = (I + (2*a2 - a1^2)*D2 + a2^2*D4);
	z = zeros(n,1);
	z(1) = 1;
	% acf of the squared process
	acf2 = AA \ z;
	%acf2 = acf2/acf2(1);
	S2   = fft(acf2);
%	S1   = sqrt(S2);
	S1   = S2;
	df   = 1/L;
	S1   = 2*S1/sum(S1*df);
end
	

