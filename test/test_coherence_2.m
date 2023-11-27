% 2023-10-02 15:40:58.741028187 +0200
% Karl Kastner, Berlin
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
%-> must be smoothed

p  = 0;
n  = 1e4;
x  = randn(n,1);
y  = p*x+(1-p)*randn(n,1);
%fx = fft(x);
%y  =fft(fy);

C = mycoherence(x,y,sqrt(n));
%plot(abs(conj(fx).*fy)/(abs(fx)*abs(fy)))

clf
plot(C)


p = 0.5;
y  = p*x+(1-p)*randn(n,1);
C = mycoherence(x,y,sqrt(n));

hold on
plot(C)

p = 1;
y  = p*x+(1-p)*randn(n,1);
C = mycoherence(x,y,sqrt(n));

hold on
plot(C)
