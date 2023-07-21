% 2023-07-18 20:54:16.706001541 +0200
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
% function a = advection_kernel(dx,dt,n,v)
function [a,out] = advection_kernel(dx,dt,n,v)
	if (length(n) == 1)
		a = zeros(n,1);
	else
		a = zeros(n);
	end

	% width of kernel
	s = dt*v./dx;

	l = floor(s);
	r = l + 1;
	p = s-l;
	% wrap integer shifts
	l = mod(l,n(1));
	r = mod(r,n(1));
	
	if (length(n)==1)
		a(l+1) = 1-p;
		a(r+1) = p;
	else
		a(l(1)+1,l(2)+1) = (1-p(1))*(1-p(2));
		a(r(1)+1,l(2)+1) = p(1)*(1-p(2));
		a(l(1)+1,r(2)+1) = (1-p(1))*p(2);
		a(r(1)+1,r(2)+1) = p(1)*p(2);
	end

	% estimate numerical diffusion coefficient
	if (nargout()>1)
		sd = sqrt(p.*(1-p))*dx;
		en = heat_equation_fundamental_std_to_time(sd,dt);		
		out.p  = p;
		out.sd = sd;
		out.en = en;	
	end
end

