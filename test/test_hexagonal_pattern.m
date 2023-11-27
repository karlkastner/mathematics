% Mon 21 Jun 15:19:21 CEST 2021
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
 
n=1e3;
f = 10;

z = generate_isotropic_pattern(f,n);

%z = z.^2;
%z = 1+(1-z).^2;
figure(2);
clf
g = graythresh(z);
g = 0.8;
subplot(2,2,1)
 imagesc(z)
 axis equal
subplot(2,2,2)
 imagesc(z>g)
 axis equal
subplot(2,2,3)
 imagesc(1-(z<1-g))
 axis equal
subplot(2,2,4);
z = z-mean(z(:));
w=1;
w = tukeywin(length(z),1);
f2 = abs(fft2(w.*z)).^2;
f2 = fftshift(f2);
f2 =f2/max(f2(:));
imagesc(log10(f2))
caxis([-4,0])
figure(3);
clf

f2 = abs(fft2(w.*z)).^2;
f2 = fftshift(f2);
%f2 = meanfilt1(f2,5);
%f2 = meanfilt1(f2',5)';
mu = r_spectrum(f2);
tmu = t_spectrum(f2);
subplot(2,2,1)
semilogy(mu/max(mu));
hold on;
subplot(2,2,2)
semilogy(tmu/max(tmu))
hold on

w = tukeywin(length(z),1);
w = w.*w';
w = w/mean(w(:));
f2 = abs(fft2(w.*z)).^2;
f2 = fftshift(f2);
mu = r_spectrum(f2);
tmu = t_spectrum(f2);
subplot(2,2,1)
semilogy(mu/max(mu))
subplot(2,2,2)
semilogy(tmu/max(tmu))

w = rwin(f2);
f2 = abs(fft2(w.*z)).^2;
f2 = fftshift(f2);
mu = r_spectrum(f2);

tmu = t_spectrum(f2);
subplot(2,2,1)
semilogy(mu/max(mu))
subplot(2,2,2)
semilogy(tmu/max(tmu))
%imagesc(w)

%semilogy(mu/max(mu))
%subplot(2,2,3)
%plot(mu)
%imagesc(w)

clf
subplot(2,2,1)
surface(10*b0,'edgecolor','none')
axis equal
colorbar
subplot(2,2,2)
plot(b0(round(end/2),:))
hold on
plot(b0(:,round(end/4)))
pause

