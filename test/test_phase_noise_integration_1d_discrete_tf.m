% Fri 27 Oct 10:45:26 CEST 2023
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
a1=0.1;
		 a2=0.1;
		 L=20;
		 n=10*L^2;
		 tf = phase_noise_integration_1d_discrete_tf(n,L,a1,a2);
		 tf_ = pni_tf(f,a1,a2);
		 e = rand(n,1);
		 b = ifft(tf.*fft(e));
		 ir=phase_noise_integration_1d_discrete_ir(n,L,a1,a2);
		 
		 subplot(2,2,1);
		 plot(ir);
		 subplot(2,2,2);
		 plot(b)

		 subplot(2,2,3);
		 cla
		 f = fourier_axis(L,n);
		 plot(f,[real(tf),real(tf_)])
%		 hold on
%		 plot(f,[real(ir),imag(tf)])
		 subplot(2,2,4)
		 plot(f,[imag(tf),imag(tf_)])

	L=40; n=2*L^2; a1=0.05; a2=0.2; x=linspace(-L/2,L/2,n); ir=fftshift(phase_noise_integration_1d_discrete_ir(n,L,a1,a2)); ir(:,2) = pni_ir(x,a1,a2); plot(x,ir)

