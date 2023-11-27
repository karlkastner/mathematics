% Mon 14 Aug 16:10:33 CEST 2023
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
function [ex,ey] = random_field_scaled(L,n,s,rs,rho,order)
	if (isempty(rs))
		rs = 1;
	end
	dx = L(1)/n(1);
	dy = L(2)/n(2);

	cx = randn(n);
	cy = randn(n);
	% enforce isotropy with rs = 1
	cy = rs*cx + sqrt(1-rs^2)*cy;
%	cx = lowpass2d_fft(x,rho,[],order);
%	cy = lowpass2d_fft(y,rho,[],order);
	% shift left neighbours
	ex = [cx(end,:) - cx(2,:)
              cx(1:end-2,:) - cx(3:end,:);
	      cx(end-1,:) - cx(1,:);
	     ];
	ey = [cy(:,end) - cy(:,2), ...
              cy(:,1:end-2) - cy(:,3:end), ...
              cy(:,end-1) - cy(:,1)];
	% n.b. : here another process, like a brownian surface could be used,
	%        though it is not straight forward how to envforce the delta constraint
	ex = s*lowpass2d_implicit(ex,rho,[],order);
	ey = s*lowpass2d_implicit(ey,rho,[],order);
	
%	ey = ex;
	
%	x0 = x;
%	y0 = y;

%	x = repmat(linspace(0,L(1),n(1))',1,n(2));
%	y = repmat(linspace(0,L(2),n(2)),n(1),1);
%	x = x + dx;
%	y = y + dy;
end
