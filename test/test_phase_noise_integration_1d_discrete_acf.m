% 2023-10-25 11:36:27.462026617 +0200
% Karl Kastner, Berlin
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
% note that in contrast to the lowpass, the sign is the opposit in front of a2 in front of the D2 term,
% and that the D1 term is present
% due to appliaction of the filter in both directions, the sign of a1 is irrelevant
		 a1 = -0.5;
		 a2 = +0.5;
		 L=100;
		 n = L^2+1;
		 [acf,A]=phase_noise_integration_1d_discrete_acf(n,L,a1,a2);
		 acf=acf/acf(1);
		x=linspace(0,L,n)';
		 y=exp(-x/5).*cos(0.89*pi/2*x);
		 y=y+[0;
		flipud(y(2:end))];
subplot(2,2,1)
		 plot(x,[acf]) , % ,y
subplot(2,2,2)
plot(real(ifft(acf)))

