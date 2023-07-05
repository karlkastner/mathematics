% 2017-04-07 14:31:30.015413021 +0200
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
%% coefficients of the Fourier series of Q|Q| 
%%	Q|Q| = Q_a^2 y					(8.5)
%%           = | cos a + cos t | (cos a + cos t)	(8.6)
%%           = a0 + a1 cos t + ... + an cos n t		(8.7)
% in general
%%            cos a is midrange
%%            cos t is tidal variation
%% c.f Dronkers 1964, eq. 8.10
% TODO this does not seem to work for alpha < 0
function an = fourier_signed_square(alpha,n)
	if (issym(alpha))
		syms pi_
	else
		pi_ = pi;
	end
	switch (n)
	case {0}
		an = 1/pi_*( (2+cos(2*alpha))*(1/2*pi_-alpha) + 3/2*sin(2*alpha));
	%	an = an/2;
	case {1}
		an = 1/pi_*(4*cos(alpha)*(1/2*pi_-alpha) + 3*sin(alpha) + 1/3*sin(3*alpha));
		%an = 1/pi_*(4*cos(alpha*(1/2*pi_-alpha)) + 3*sin(alpha) + 1/3*sin(3*alpha));
	case {2}
		an = 1/pi_*(1/2*pi_ - alpha + 2/3*sin(2*alpha) - 1/12*sin(4*alpha));
	otherwise % note: this does not apply to 1..2, division by zero!
	an = (-1)^(n+1)*2/(n*pi_)...
		*( sin((n-2)*alpha)/((n-1)*(n-2)) ...
		   - 2*sin(n*alpha)/((n+1)*(n-1)) ...
		   + sin((n+2)*alpha)/((n+2)*(n+1)) );
	end
end

