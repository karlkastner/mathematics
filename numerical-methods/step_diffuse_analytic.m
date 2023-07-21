% Mon  2 May 14:18:38 CEST 2022
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
%% analytic solution to the heat equation
%
%% the spectral solution is not positivity preserving as it results in 
%% spurious oscillations, this is avoided here, by integrating over segments
%% rather than sampling at gridpoints
%
% dy/dt = e*D^2*y
function z = diffuse_analytic(t,z,nx,dx,e)


	% diffuse along first dimension
	z = diffuse(e(1)*t,z,dx(1));
	
	% the kernel is gaussian, so the dimensions are separable
	% convolve along x
	if (length(dx)>1)
		% convolve along y
		z = diffuse(e(2)*t,z.',dx(2)).';
	end

	%heat_equation_fundamental_solution(dt,x,d,t0,x0);
	% g = ifft(T*fft(b))
	% fft(g) = T*fft(b)

function z = diffuse_()
	% standard deviation of the initial unit impulse
	sd = sqrt(12*dx);
	t0 = heat_equation_std_to_time(sd,te);
	g  = heat_equation_fundamental_solution(te,x,1,t0,0);
	fg = fft(g)
	fz = fft(z);
	z = ifft(fg.*fz);
end
	
% n.b., this 
function z = diffuse(et,z,dx)
	%nx = size(z,1);
	%dx = L/nx;

	if (0==mod(length(z),2))
		x = dx*[0:nx/2-1,-nx/2:-1]';
	else
		x = dx*[0:(nx-3)/2,-(nx+1)/2:-1]';
	end
	if (0)
		% fw euler:
		% z = (I + dt*e*D2)^k*z
		% z = (I - dt*e*2/dx^2)*z

		% sample at midpoints
		% y = normpdf(x,x0,sqrt(2*d*(t+t0)));
		g = dx./sqrt(4*pi*et)*exp(-(x.^2)./(4*et));
		% g(0) ~ 1/
	else
		% integrated analytically
		g = dx/2*(erf((x+dx/2)/sqrt(4*et)) - erf((x-dx/2)/sqrt(4*et)));
%	g = dx./sqrt(4*pi*t)*exp(-((x-0.5*dx).^2)./(4*t));
%	g = g + dx./sqrt(4*pi*t)*exp(-((x+0.5*dx).^2)./(4*t));
%	g = g/2;
	end
	mg = max(g);
	if (mg>1)
		warning('max(g)>1');
	end
	g = g/sum(g);
	if (0)
		z = cconv(z,g,length(z));
	else
%hold on
%plot(abs(fft(g).^100))
		z = ifft(fft(g).*fft(z));
	end	
end % diffuse	
end % diffuse_analytic

