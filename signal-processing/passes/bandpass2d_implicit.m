% Thu 24 Jun 14:39:10 CEST 2021
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
%
%% bandpass filter the surface x by solving the implicit relation:
%
% function [y] = bandpass2d_implicit(x,rho,a,order,direct)
%
% y[i] = rho*(y[i-1] + y[i+1]) + (1-2*rho)*x[i]
% (1 - 2 rho)y = rho*[1,-2,1]*y
% ((1-4*rho)*I - rho*D2)*y = (1-4*rho)*x
% y = (1-4*rho)*(((1-4*rho)*I - rho*D2))^-1*x
% lowpass:
% A = (((1-2*(rhox+rhoy))*I - rhox*D2x - rhoy*D2y));
%
%
function [y] = bandpass2d_implicit(x,rho,a,order,direct)
	if (any(rho>=1))
		warning('rho must be smaller 1');
	end
	if (length(rho) < 2)
		rho = [rho,rho];
	end
	if (nargin()<5||isempty(direct))
		direct = false;
	end
	n  = size(x);
	tol = [];
	maxit = 2*max(n);
	bc = {'circular','circular'};
	[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,n-1,2,bc);
	I = speye(prod(n));

	% TODO is the normalization correct?
	r  = rho./(1 - 2*rho + rho.*rho)
	D2 = (r(1)*D2x + r(2)*D2y);

	if (nargin()>2&&~isempty(a))
		s = sin(a);
		c = cos(a);
		D2 = c*c*D2x + 2*s*c*Dxy + s*s*D2y;
		% TODO asymmetric scaling
		D2 = r(1)*D2;
	end

	A = (I - D2);
	y = flat(x);
	if (mod(order,1)==0)
	for idx=1:order
	if (direct)
		A2 = A*A;
		y = (A2 \ ((A-I)*(y)));
	else
		y =     pcg(A, y,tol,maxit);
		y = y - pcg(A, y,tol,maxit);
		if (0)
			% implementation in two steps is faster, better conditioned!
			A2 = A*A;
			y = pcg(A2,(A-I)*flat(y),tol,maxit,[],[],flat(x));
		end
	end
		%y = y/norm(y);
	end
	else
		A2 = A*A;
%		[V,E]=eigs(
	end
	y = reshape(y,n);
end

