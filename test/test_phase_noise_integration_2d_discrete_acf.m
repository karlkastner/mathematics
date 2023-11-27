% 2023-10-26 16:07:49.087053730 +0200 
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
		 a1 = -0.05;
		 a2 = +0.05;
		 a3 = -0.7;
		 L  = 20*[1,1];
		 n  = L.^2;
	%	 [acf,A]=phase_noise_integration_2d_discrete_acf(n,L,a1,a2,a3);
		 [tf,A]=phase_noise_integration_2d_discrete_tf(n,L,a1,a2,a3);
		S = abs(tf).^2;
		%x=linspace(0,L,n)';
		% y=exp(-x/5).*cos(0.89*pi/2*x);
		 %y=y+[0;
		%flipud(y(2:end))];

		rng(0)	
		e = randn(n);
		%T = sqrt(S);
		b = real(ifft2(tf.*fft2(e)));


%figure(1)
%subplot(2,3,1)
%		imagesc(acf)
%subplot(2,3,2)
%	plot(acf(1,:));
%subplot(2,3,3)
%	plot(acf(:,1));

		 %plot(x,[acf]) , % ,y
%S = real(ifft2(acf));

%subplot(2,3,4)
%imagesc(fftshift(log(S)))
%imagesc(fftshift((S)))

figure(2)
subplot(2,3,1)
imagesc(real(b))
axis equal
axis tight

subplot(2,3,2)
bt = quantile(flat(b),0.66);
imagesc(real(b)>bt)
axis equal
axis tight

subplot(2,3,3)
[fx,fy] = fourier_axis_2d(n,L);
imagesc(fftshift(fx),fftshift(fy),fftshift(abs(tf).^2))
axis equal
axis tight

subplot(2,3,4)
plot(mean(S))

subplot(2,3,5)
plot(mean(S,2))


