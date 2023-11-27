% 2021-06-16 14:17:30.331659100 +0200
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

%A = [   0,    -x(1); x(2),  0]
%eig(A)
%ydot = @(t,x) [0, -x(1); x(2), 0]*x + [0.1; 0]
%ydot = @(t,x) [-0.05*x(1)*x(2); 0.0125/4*x(1)*x(2)] + [4*1.9*x(1); -2*x(2)]/2
%[0.1*sign(x(1)-1).*(x(1)-1); 0];
%ydot = @(t,x) [0,-1; 1, 0]*x + (1-hypot(x(1),x(2))).*x


if (~exist('y','var'))
y0 = [ 1, 0];
T  = [ 0,2*pi*100];
clf
%max(y)
[t,y] = nonlinear_oscillator_noisy(T,y0,0,0);
%plot(t,y);
end

if (1)

clf
T = linspace(0,T(2),100*length(t));
[t,y] = oscillator_noisy(T,y0,1,0);
t=t/(2*pi);
hold on

clf
subplot(2,1,1)
plot(t,y(:,1))
hold on
subplot(2,2,3)
plot(t,hypot(y(:,1),y(:,2)))
hold on
subplot(2,2,4)
%plot(t,
p = unwrap(atan2(y(:,2),y(:,1)));
%);
%hold on
y_ = y;

s = 10*1e-2/sqrt(t(2)-t(1));
c=1;
%T = linspace(0,T(end),length(t));
[t,y] = oscillator_noisy(T,y0,1,0);
%[t,y] = euler(ydot,T,y0);
t=t/(2*pi);
subplot(2,1,1)
plot(t,y(:,1))
subplot(2,2,3)
h=hypot(y(:,1),y(:,2));
plot(t,h);
ylim([0,1.1])
subplot(2,2,4)
p(:,2) = unwrap(atan2(y(:,2),y(:,1)));
dp=p(:,2)-p(:,1);
plot(t,[dp,[diff(dp);0]])
hold on
rms(h-1)

figure(3);
clf

y=[y_(:,1),y(:,1),y_(:,1)+s*randn(size(y_,1),1)];
y = y-mean(y);
F = fft(y);
F = abs(F).^2;
F = F./norm(F);
subplot(2,1,1)
loglog(F)
subplot(2,1,2)
loglog(trifilt1(F,11))


end
