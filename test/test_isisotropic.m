% Fri 10 Nov 13:06:13 CET 2023
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
>> for idx=1:10; n=1000; m=50*(idx-1)+1; t = 2*pi*(0:n-1)'/n; S = zeros(n,1); S(m) = 1; S(n/2+m) = 1; S_(:,1)=1; [piso,t0,c,S_] = isotropy(t,S); piso, t0, subplot(2,5,idx); plot(t,[S,0.5+0.5*(S_-mean(S_))./max(S_-mean(S_)),0.5+0.5*sin(2*t+t0)]),c(2:3)./hypot(c(2),c(3)); ylim([0,1]);end
n=600; t=(0:n-1)/n*2*pi; S = zeros(n,1); S(1:n/6:n-1) = 1; isisotropic(t,S); plot(S)
n=600; t=(0:n-1)/n*2*pi; S = zeros(n,1); S(1:n/6:n-1) = 1/3; S(1)=1; S(n/2+1)=1; isisotropic(t,S), plot(S)
